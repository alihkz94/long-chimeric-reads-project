"""
This script processes BLAST XML files and corresponding FASTA files that they are chimeric (multi_alignment), 
summarizing their analysis results and producing both detailed outcomes and summary statistics.
"""

import os
import glob
import multiprocessing as mp
from Bio import SeqIO
from Bio.Blast import NCBIXML
import pandas as pd
import logging
from dataclasses import dataclass
from typing import List, Dict, Tuple, Generator
from collections import defaultdict
import itertools
import gc
import numpy as np
from contextlib import contextmanager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

@dataclass
class AlignmentSegment:
    """Represents a single alignment segment from BLAST results"""
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
    """Context manager to ensure proper resource cleanup"""
    try:
        yield
    finally:
        gc.collect()

class ChimeraAnalyzer:
    def __init__(self, xml_file: str, fasta_file: str, chunk_size: int = 1000):
        self.xml_file = xml_file
        self.fasta_file = fasta_file
        self.chunk_size = chunk_size
        
    def _load_sequences_generator(self) -> Generator[Dict[str, str], None, None]:
        """Generator to load sequences in chunks to conserve memory"""
        current_chunk = {}
        for i, record in enumerate(SeqIO.parse(self.fasta_file, "fasta"), 1):
            current_chunk[record.id] = str(record.seq)
            if i % self.chunk_size == 0:
                yield current_chunk
                current_chunk = {}
        if current_chunk:
            yield current_chunk

    def _extract_taxonomy(self, hit_def: str) -> str:
        try:
            return hit_def.split(';', 1)[1].strip()
        except IndexError:
            return "Unknown taxonomy"

    def _analyze_segment_relationships(self, segments: List[AlignmentSegment], query_id: str, 
                                    query_length: int) -> List[Dict]:
        """Analyze relationships between all segments"""
        segment_data = []
        
        if not segments:
            return segment_data

        # Sort segments by start position
        segments.sort(key=lambda x: x.start)
        total_coverage = sum(seg.coverage for seg in segments)
        
        # Generate all possible adjacent pairs and triplets
        for i in range(len(segments) - 1):
            current_seg = segments[i]
            next_seg = segments[i + 1]
            
            # Calculate overlap or gap between segments
            overlap = current_seg.end - next_seg.start
            gap = next_seg.start - current_seg.end if overlap < 0 else 0
            
            segment_info = {
                # Basic sequence information
                'sequence_id': query_id,
                'sequence_length': query_length,
                'total_segments': len(segments),
                'segment_pair': f"{i+1}-{i+2}",
                
                # Breakpoint information
                'breakpoint_position': (current_seg.end + next_seg.start) // 2,
                'overlap_size': max(0, overlap),
                'gap_size': gap,
                
                # First segment characteristics
                'segment1_start': current_seg.start,
                'segment1_end': current_seg.end,
                'segment1_length': current_seg.end - current_seg.start + 1,
                'segment1_coverage': current_seg.coverage,
                'segment1_identity': current_seg.identity,
                'segment1_taxonomy': current_seg.taxonomy,
                'segment1_strand': current_seg.strand,
                'segment1_evalue': current_seg.e_value,
                
                # Second segment characteristics
                'segment2_start': next_seg.start,
                'segment2_end': next_seg.end,
                'segment2_length': next_seg.end - next_seg.start + 1,
                'segment2_coverage': next_seg.coverage,
                'segment2_identity': next_seg.identity,
                'segment2_taxonomy': next_seg.taxonomy,
                'segment2_strand': next_seg.strand,
                'segment2_evalue': next_seg.e_value,
                
                # Overall characteristics
                'total_coverage': total_coverage,
                'different_taxa': current_seg.taxonomy != next_seg.taxonomy,
                'different_strands': current_seg.strand != next_seg.strand
            }
            
            # If there's a third segment, include its information
            if i + 2 < len(segments):
                third_seg = segments[i + 2]
                segment_info.update({
                    'segment3_start': third_seg.start,
                    'segment3_end': third_seg.end,
                    'segment3_length': third_seg.end - third_seg.start + 1,
                    'segment3_coverage': third_seg.coverage,
                    'segment3_identity': third_seg.identity,
                    'segment3_taxonomy': third_seg.taxonomy,
                    'segment3_strand': third_seg.strand,
                    'segment3_evalue': third_seg.e_value
                })
            
            segment_data.append(segment_info)
            
        return segment_data

    def process_chunk(self, xml_chunk: List[str]) -> List[Dict]:
        """Process a chunk of BLAST records"""
        chunk_data = []
        
        for record in xml_chunk:
            query_id = record.query
            query_length = record.query_length
            
            if not record.alignments:
                continue
            
            segments = []
            for alignment in record.alignments:
                hit_id = alignment.hit_def
                taxonomy = self._extract_taxonomy(hit_id)
                
                for hsp in alignment.hsps:
                    coverage = ((hsp.query_end - hsp.query_start + 1) / query_length) * 100
                    identity = (hsp.identities / hsp.align_length) * 100
                    strand = 1 if hsp.strand == ('Plus', 'Plus') else -1
                    
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
            
            if segments:
                chunk_data.extend(self._analyze_segment_relationships(segments, query_id, query_length))
        
        return chunk_data

def process_file(args: Tuple[str, str, int]) -> pd.DataFrame:
    """Process a single file pair (for multiprocessing)"""
    xml_file, fasta_file, chunk_size = args
    
    try:
        analyzer = ChimeraAnalyzer(xml_file, fasta_file, chunk_size)
        all_data = []
        
        with open(xml_file) as result_handle:
            # Process BLAST records in chunks
            records = NCBIXML.parse(result_handle)
            while True:
                chunk = list(itertools.islice(records, chunk_size))
                if not chunk:
                    break
                
                with managed_resource():
                    chunk_data = analyzer.process_chunk(chunk)
                    if chunk_data:
                        all_data.extend(chunk_data)
        
        if all_data:
            df = pd.DataFrame(all_data)
            df['source_xml'] = os.path.basename(xml_file)
            df['source_fasta'] = os.path.basename(fasta_file)
            return df
        
    except Exception as e:
        logging.error(f"Error processing {xml_file}: {e}")
    
    return pd.DataFrame()

def process_multiple_files(xml_files: List[str], fasta_files: List[str], 
                         output_dir: str, chunk_size: int = 1000):
    """Process multiple BLAST XML and FASTA file pairs using multiprocessing"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Determine number of processes
    num_processes = max(1, mp.cpu_count()) #can put "-1" after mp.cpu_count() or any other numbers to decrese the number of processes
    logging.info(f"Using {num_processes} processes")
    
    # Prepare arguments for multiprocessing
    args_list = [(xml, fasta, chunk_size) for xml, fasta in zip(xml_files, fasta_files)]
    
    # Process files in parallel
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(process_file, args_list)
    
    # Combine results
    final_df = pd.concat(results, ignore_index=True)
    
    if not final_df.empty:
        # Save detailed results in chunks
        chunk_size = 50000  # Adjust based on your memory constraints
        for i, chunk_df in enumerate(np.array_split(final_df, max(1, len(final_df) // chunk_size))):
            output_file = os.path.join(output_dir, f"chimera_analysis_part{i+1}.csv")
            chunk_df.to_csv(output_file, index=False)
        
        # Generate summary statistics
        summary_stats = pd.DataFrame({
            'total_sequences': [len(final_df['sequence_id'].unique())],
            'total_breakpoints': [len(final_df)],
            'max_segments_per_sequence': [final_df['total_segments'].max()],
            'avg_sequence_length': [final_df['sequence_length'].mean()],
            'avg_segment_coverage': [final_df[['segment1_coverage', 'segment2_coverage']].mean().mean()],
            'avg_segment_identity': [final_df[['segment1_identity', 'segment2_identity']].mean().mean()],
            'different_taxa_count': [final_df['different_taxa'].sum()],
            'different_strands_count': [final_df['different_strands'].sum()],
            'multi_segment_sequences': [(final_df['total_segments'] > 2).sum()]
        })
        
        summary_file = os.path.join(output_dir, "analysis_summary.csv")
        summary_stats.to_csv(summary_file, index=False)
        
        logging.info(f"Analysis complete. Results saved to {output_dir}")
        return final_df
    
    return None

def get_matching_files(xml_dir: str, fasta_dir: str) -> Tuple[List[str], List[str]]:
    """
    Find matching FASTA and XML files based on their ERR numbers.
    Returns two lists of file paths (xml_files, fasta_files) where indexes correspond to matching pairs.
    """
    # Get all files
    fasta_files = glob.glob(os.path.join(fasta_dir, "*.fasta"))
    xml_files = glob.glob(os.path.join(xml_dir, "*.xml"))
    
    # Extract ERR numbers and create mapping
    fasta_map = {os.path.basename(f).split('_')[0]: f for f in fasta_files}
    xml_map = {os.path.basename(f).split('_')[0]: f for f in xml_files}
    
    # Find common ERR numbers
    common_errs = sorted(set(fasta_map.keys()) & set(xml_map.keys()))
    
    # Create matched lists
    matched_xml_files = [xml_map[err] for err in common_errs]
    matched_fasta_files = [fasta_map[err] for err in common_errs]
    
    return matched_xml_files, matched_fasta_files

def main():
    # Directories containing XML and FASTA files
    xml_dir = "../xml_dada2_nonchimeras"  # Replace with your XML files directory
    fasta_dir = "."  # Replace with your FASTA files directory
    output_dir = "chimera_analysis_output"

    # Get matching files
    xml_files, fasta_files = get_matching_files(xml_dir, fasta_dir)

    # Log the matched files
    logging.info(f"Found {len(xml_files)} matching file pairs")
    for xml, fasta in zip(xml_files, fasta_files):
        logging.info(f"Processing pair: {os.path.basename(xml)} - {os.path.basename(fasta)}")

    # Process the files
    results = process_multiple_files(xml_files, fasta_files, output_dir)

    if results is not None:
        logging.info("Analysis completed successfully")
    else:
        logging.error("No results were generated")

if __name__ == "__main__":
    main()
