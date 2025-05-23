#!/usr/bin/env python3
"""
otu_taxonomy_extractor.py

A high-performance, multiprocessing, streaming BLAST-XML taxonomy & hit-info extractor
that matches sequences from OTUs.fasta with sequences in BLAST XML files
to extract taxonomic information and key HSP metrics.

Requirements:
    - Python 3.6+
    - lxml       (pip install lxml)

Usage:
    python otu_taxonomy_extractor.py \
        --input . \
        --fasta OTUs.fasta \
        --output taxonomy.tsv \
        --workers 8
"""

import os
import argparse
import logging
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from lxml import etree

def normalize_sequence(seq):
    """Normalize a sequence for comparison (uppercase, no whitespace)."""
    return re.sub(r'\s+', '', seq.upper())

def parse_fasta(fasta_path):
    """
    Parse FASTA file and return two dicts:
      id_to_seq: {sequence_id: sequence}
      seq_to_id: {normalized_sequence: sequence_id}
    """
    id_to_seq = {}
    seq_to_id = {}
    with open(fasta_path) as f:
        curr_id = None
        curr_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if curr_id:
                    seq = ''.join(curr_seq)
                    id_to_seq[curr_id] = seq
                    seq_to_id[normalize_sequence(seq)] = curr_id
                curr_id = line[1:].split()[0]
                curr_seq = []
            elif curr_id:
                curr_seq.append(line)
        # last
        if curr_id:
            seq = ''.join(curr_seq)
            id_to_seq[curr_id] = seq
            seq_to_id[normalize_sequence(seq)] = curr_id
    return id_to_seq, seq_to_id

def is_self_hit(hit_def: str) -> bool:
    """Detect BLAST self-hit by ';size=' in the def."""
    return ';size=' in (hit_def or '')

def extract_sequences_from_xml(xml_path):
    """
    First pass: stream-parse XML and extract, per query sequence,
    the first non-self hit’s taxonomy + HSP fields.
    Returns { norm_query_seq: { all fields… } }
    """
    seq_tax_map = {}
    for _, it in etree.iterparse(xml_path, events=('end',), tag='Iteration'):
        try:
            for hit in it.findall('Iteration_hits/Hit'):
                hdef = hit.findtext('Hit_def')
                if not hdef or is_self_hit(hdef):
                    continue
                # pull hit-level fields
                hit_num       = hit.findtext('Hit_num')
                hit_id_text   = hit.findtext('Hit_id')
                hit_acc       = hit.findtext('Hit_accession')
                hit_len       = hit.findtext('Hit_len')
                # first Hsp
                hsp = hit.find('Hit_hsps/Hsp')
                if hsp is None:
                    continue
                qseq = hsp.findtext('Hsp_qseq')
                if not qseq:
                    continue
                norm = normalize_sequence(qseq)
                # extract all Hsp fields
                bit_score = float(hsp.findtext('Hsp_bit-score') or 0)
                score     = float(hsp.findtext('Hsp_score') or 0)
                evalue    = float(hsp.findtext('Hsp_evalue') or 0)
                qfrom     = int(hsp.findtext('Hsp_query-from') or 0)
                qto       = int(hsp.findtext('Hsp_query-to') or 0)
                hfrom     = int(hsp.findtext('Hsp_hit-from') or 0)
                hto       = int(hsp.findtext('Hsp_hit-to') or 0)
                qframe    = int(hsp.findtext('Hsp_query-frame') or 0)
                hframe    = int(hsp.findtext('Hsp_hit-frame') or 0)
                identity  = int(hsp.findtext('Hsp_identity') or 0)
                positive  = int(hsp.findtext('Hsp_positive') or 0)
                gaps      = int(hsp.findtext('Hsp_gaps') or 0)
                align_len = int(hsp.findtext('Hsp_align-len') or 0)

                seq_tax_map[norm] = {
                    'hit_num':       hit_num,
                    'hit_id':        hit_id_text,
                    'hit_def':       hdef,
                    'hit_accession': hit_acc,
                    'hit_len':       hit_len,
                    'hsp_num':       hsp.findtext('Hsp_num'),
                    'bit_score':     round(bit_score, 2),
                    'score':         round(score, 2),
                    'evalue':        round(evalue, 2),
                    'qfrom':         qfrom,
                    'qto':           qto,
                    'hfrom':         hfrom,
                    'hto':           hto,
                    'qframe':        qframe,
                    'hframe':        hframe,
                    'identity':      identity,
                    'positive':      positive,
                    'gaps':          gaps,
                    'align_len':     align_len,
                    # coverage & identity%
                    'coverage':      None,  # filled later
                    'identity_pct':  None
                }
                break
        except Exception as e:
            logging.warning(f"Error parsing {xml_path}: {e}")
        finally:
            it.clear()
    return seq_tax_map

def parse_xml_stream(xml_path, otu_seq_to_id, id_to_seq):
    """
    Match each extracted seq → hit_info dict to OTU IDs.
    Compute coverage & identity% from align_len and original OTU length.
    Returns { otu_id: hit_info_dict }
    """
    mapping = {}
    seq_tax = extract_sequences_from_xml(xml_path)
    for norm_seq, info in seq_tax.items():
        otu_id = otu_seq_to_id.get(norm_seq)
        if not otu_id:
            continue
        length = len(id_to_seq[otu_id])
        cov = (info['align_len'] / length * 100) if length else 0
        info['coverage']     = round(cov, 2)
        info['identity_pct'] = round((info['identity'] / info['align_len'] * 100) if info['align_len'] else 0, 2)
        mapping[otu_id] = info
    logging.info(f"{xml_path}: matched {len(mapping)} OTUs")
    return mapping

def gather_xmls(path):
    if os.path.isdir(path):
        return sorted(
            os.path.join(path, f) for f in os.listdir(path)
            if f.lower().endswith(('.xml','.blastxml'))
        )
    elif os.path.isfile(path):
        return [path]
    else:
        raise FileNotFoundError(path)

def main():
    p = argparse.ArgumentParser(description="OTU-based BLAST XML extractor")
    p.add_argument('-i','--input',  required=True, help="XML file or dir")
    p.add_argument('-f','--fasta',  required=True, help="OTUs.fasta")
    p.add_argument('-o','--output', required=True, help="Output TSV (tab-delimited)")
    p.add_argument('-w','--workers',type=int,default=8,
                   help="Number of processes (default=8)")
    p.add_argument('-v','--verbose',action='store_true')
    args = p.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s"
    )

    # FASTA → dicts
    logging.info(f"Reading FASTA {args.fasta}")
    id_to_seq, seq_to_id = parse_fasta(args.fasta)

    # XML files
    xmls = gather_xmls(args.input)
    logging.info(f"Found {len(xmls)} XML(s)")

    # Parallel extraction
    combined = {}
    with ProcessPoolExecutor(max_workers=args.workers) as exe:
        futs = {
            exe.submit(parse_xml_stream, xf, seq_to_id, id_to_seq): xf
            for xf in xmls
        }
        for fut in as_completed(futs):
            xf = futs[fut]
            try:
                combined.update(fut.result())
                logging.info(f"✅ {xf}")
            except Exception as e:
                logging.error(f"❌ {xf}: {e}")

    # Write out tab-delimited table with only essential fields
    fields = ['hit_def', 'evalue', 'identity_pct', 'coverage']
    with open(args.output, 'w') as out:
        out.write('OTU_ID\t' + '\t'.join(fields) + '\n')
        for otu, info in sorted(combined.items()):
            vals = [str(info[f]) for f in fields]
            out.write(otu + '\t' + '\t'.join(vals) + '\n')

    logging.info(f"Wrote {len(combined)} records to {args.output}")

if __name__ == '__main__':
    main()
