"""
This script processes BLAST XML files and corresponding FASTA files that they are chimeric (multi_alignment), 
summarizing their analysis results and producing both detailed outcomes and summary statistics.
"""

import os
import glob
import multiprocessing as mp
from Bio.Blast import NCBIXML
import pandas as pd
import logging
from dataclasses import dataclass
from typing import List, Dict, Tuple
import itertools
import gc
from contextlib import contextmanager

# Configure logging for tracking progress and errors
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

@dataclass
class AlignmentSegment:
    """
    Represents a single alignment segment (HSP) from BLAST results.
    Each instance includes positional, identity, coverage, taxonomic, and other metadata.
    """
    start: int
    end: int
    identity: float
    coverage: float
    hit_id: str
    taxonomy: str
    score: float
    e_value: float
    strand: int

@contextmanager
def managed_resource():
    """
    Context manager for resource management to handle garbage collection explicitly.
    Ensures memory cleanup after processing.
    """
    try:
        yield
    finally:
        gc.collect()

class ChimeraAnalyzer:
    """
    Analyzes BLAST XML results and associates hits with sample or database taxonomy.
    Provides functionality for identifying and processing alignment segments.
    """
    def __init__(self, xml_file: str, fasta_file: str, chunk_size: int = 1000):
        self.xml_file = xml_file
        self.fasta_file = fasta_file
        self.chunk_size = chunk_size

    def _get_query_id_from_header(self, header: str) -> str:
        """
        Extracts the query ID (base part) from the FASTA header.
        This is used for identifying self-hits and sample-specific hits.
        """
        return header.split(';')[0]

    def _is_self_hit(self, query_id: str, hit_def: str) -> bool:
        """
        Determines if the BLAST hit is a self-hit by comparing the query ID and hit definition.
        Self-hits are ignored in downstream processing.
        """
        hit_id = hit_def.split(';')[0]
        query_base_id = self._get_query_id_from_header(query_id)
        return hit_id == query_base_id

    def _is_database_hit(self, hit_def: str) -> bool:
        """
        Identifies if the hit belongs to a reference database based on specific taxonomic prefixes.
        Database hits typically include standardized prefixes (e.g., k__, p__, c__).
        """
        return any(prefix in hit_def for prefix in ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'])

    def _extract_taxonomy(self, hit_def: str, is_sample_hit: bool = False) -> str:
        """
        Extracts taxonomic information from the hit definition.
        For sample hits, the entire hit definition is returned; for database hits, only taxonomy is extracted.
        """
        if is_sample_hit:
            return hit_def
        try:
            tax_part = hit_def.split(';', 1)[1].strip()
            return tax_part if tax_part else "Unknown taxonomy"
        except IndexError:
            return "Unknown taxonomy"

    def process_chunk(self, xml_chunk: List[str]) -> List[Dict]:
        """
        Processes a chunk of BLAST XML records and extracts segment information.
        Identifies sample and database hits and calculates coverage and identity for each HSP.
        """
        chunk_data = []
        for record in xml_chunk:
            query_id = record.query
            query_length = record.query_length
            
            if not record.alignments:  # Skip records with no alignments
                continue
            
            # Separate hits into database and sample categories
            database_hits = []
            sample_hits = []
            for alignment in record.alignments:
                if self._is_self_hit(query_id, alignment.hit_def):  # Ignore self-hits
                    continue
                if self._is_database_hit(alignment.hit_def):
                    database_hits.append(alignment)
                else:
                    sample_hits.append(alignment)
            
            # Prioritize database hits, fall back to sample hits if none exist
            hits_to_process = database_hits if database_hits else sample_hits
            segments = []
            if hits_to_process:
                alignment = hits_to_process[0]  # Process the first alignment
                is_sample_hit = alignment in sample_hits
                hit_id = alignment.hit_def
                taxonomy = self._extract_taxonomy(hit_id, is_sample_hit)
                
                # Process HSPs for this alignment
                for hsp in alignment.hsps:
                    coverage = ((hsp.query_end - hsp.query_start + 1) / query_length) * 100
                    identity = (hsp.identities / hsp.align_length) * 100
                    strand = 1 if hsp.strand == ('Plus', 'Plus') else -1
                    
                    # Create an AlignmentSegment for each HSP
                    segment = AlignmentSegment(
                        start=hsp.query_start,
                        end=hsp.query_end,
                        identity=identity,
                        coverage=coverage,
                        hit_id=hit_id,
                        taxonomy=taxonomy,
                        score=hsp.score,
                        e_value=hsp.expect,
                        strand=strand
                    )
                    segments.append(segment)
            
            # Analyze relationships between identified segments
            if segments:
                chunk_data.extend(self._analyze_segment_relationships(segments, query_id, query_length))
        
        return chunk_data

    def _analyze_segment_relationships(self, segments: List[AlignmentSegment], query_id: str, 
                                    query_length: int) -> List[Dict]:
        """
        Analyzes relationships between adjacent segments in terms of overlap, gaps, and taxonomic differences.
        Provides detailed information for potential chimera analysis.
        """
        segment_data = []
        if not segments:
            return segment_data

        # Sort segments by start position for structured analysis
        segments.sort(key=lambda x: x.start)
        total_coverage = sum(seg.coverage for seg in segments)
        
        # Iterate through adjacent segment pairs
        for i in range(len(segments) - 1):
            current_seg = segments[i]
            next_seg = segments[i + 1]
            
            # Calculate overlap or gap between segments
            overlap = current_seg.end - next_seg.start
            gap = next_seg.start - current_seg.end if overlap < 0 else 0
            
            # Collect detailed information about segment pairs
            segment_info = {
                'sequence_id': query_id,
                'sequence_length': query_length,
                'total_segments': len(segments),
                'segment_pair': f"{i+1}-{i+2}",
                'breakpoint_position': (current_seg.end + next_seg.start) // 2,
                'overlap_size': max(0, overlap),
                'gap_size': gap,
                'segment1_start': current_seg.start,
                'segment1_end': current_seg.end,
                'segment1_length': current_seg.end - current_seg.start + 1,
                'segment1_coverage': current_seg.coverage,
                'segment1_identity': current_seg.identity,
                'segment1_taxonomy': current_seg.taxonomy,
                'segment1_strand': current_seg.strand,
                'segment1_evalue': current_seg.e_value,
                'segment2_start': next_seg.start,
                'segment2_end': next_seg.end,
                'segment2_length': next_seg.end - next_seg.start + 1,
                'segment2_coverage': next_seg.coverage,
                'segment2_identity': next_seg.identity,
                'segment2_taxonomy': next_seg.taxonomy,
                'segment2_strand': next_seg.strand,
                'segment2_evalue': next_seg.e_value,
                'total_coverage': total_coverage,
                'different_taxa': current_seg.taxonomy != next_seg.taxonomy,
                'different_strands': current_seg.strand != next_seg.strand
            }
            segment_data.append(segment_info)
        
        return segment_data

def process_file(args: Tuple[str, str, str, int]) -> None:
    """Process a single file pair and save results directly"""
    xml_file, fasta_file, output_dir, chunk_size = args
    
    try:
        # Extract base name for output file
        base_name = os.path.basename(xml_file).split('_')[0]
        output_file = os.path.join(output_dir, f"{base_name}_sequence_details.csv")
        
        analyzer = ChimeraAnalyzer(xml_file, fasta_file, chunk_size)
        
        # Process in chunks and write directly to file
        with open(xml_file) as result_handle:
            records = NCBIXML.parse(result_handle)
            first_chunk = True
            
            while True:
                chunk = list(itertools.islice(records, chunk_size))
                if not chunk:
                    break
                
                with managed_resource():
                    chunk_data = analyzer.process_chunk(chunk)
                    if chunk_data:
                        df = pd.DataFrame(chunk_data)
                        df['source_xml'] = os.path.basename(xml_file)
                        df['source_fasta'] = os.path.basename(fasta_file)
                        
                        # Write mode: 'w' for first chunk, 'a' for subsequent chunks
                        mode = 'w' if first_chunk else 'a'
                        header = first_chunk
                        df.to_csv(output_file, mode=mode, header=header, index=False)
                        first_chunk = False
        
        logging.info(f"Completed processing {xml_file}")
        
    except Exception as e:
        logging.error(f"Error processing {xml_file}: {e}")

def process_multiple_files(xml_files: List[str], fasta_files: List[str], 
                         output_dir: str, chunk_size: int = 1000):
    """Process multiple BLAST XML and FASTA file pairs using multiprocessing"""
    os.makedirs(output_dir, exist_ok=True)
    
    num_processes = max(1, mp.cpu_count()) # Use all available cores, if you want to modify put -1 or -2 before the last paranteses
    logging.info(f"Using {num_processes} processes")
    
    args_list = [(xml, fasta, output_dir, chunk_size) 
                 for xml, fasta in zip(xml_files, fasta_files)]
    
    with mp.Pool(processes=num_processes) as pool:
        pool.map(process_file, args_list)
    
    logging.info(f"Analysis complete. Results saved to {output_dir}")

def get_matching_files(xml_dir: str, fasta_dir: str) -> Tuple[List[str], List[str]]:
    """Find matching FASTA and XML files based on their ERR numbers"""
    fasta_files = glob.glob(os.path.join(fasta_dir, "*.fasta"))
    xml_files = glob.glob(os.path.join(xml_dir, "*.xml"))
    
    fasta_map = {os.path.basename(f).split('_')[0]: f for f in fasta_files}
    xml_map = {os.path.basename(f).split('_')[0]: f for f in xml_files}
    
    common_errs = sorted(set(fasta_map.keys()) & set(xml_map.keys()))
    
    matched_xml_files = [xml_map[err] for err in common_errs]
    matched_fasta_files = [fasta_map[err] for err in common_errs]
    
    return matched_xml_files, matched_fasta_files

def main():
    xml_dir = "."
    fasta_dir = "./multi"
    output_dir = "chimera_analysis_output"
    
    xml_files, fasta_files = get_matching_files(xml_dir, fasta_dir)
    
    logging.info(f"Found {len(xml_files)} matching file pairs")
    for xml, fasta in zip(xml_files, fasta_files):
        logging.info(f"Processing pair: {os.path.basename(xml)} - {os.path.basename(fasta)}")
    
    process_multiple_files(xml_files, fasta_files, output_dir)

if __name__ == "__main__":
    main()
