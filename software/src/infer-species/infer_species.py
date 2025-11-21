#!/usr/bin/env python3
"""
Species Inference Script for Single-Cell RNA-seq Data

This script analyzes count matrices from single-cell RNA-seq data to infer
the species based on gene identifiers. It supports various gene identifier
formats including gene symbols, Ensembl IDs, and Entrez IDs.

Usage:
    python infer_species.py <input_file>
    
Output:
    Creates a species.txt file containing the inferred species name.
"""

import argparse
import sys
import re
import random
from pathlib import Path
from typing import List, Tuple, Dict
import polars as pl
import scanpy as sc


class SpeciesInferenceError(Exception):
    """Custom exception for species inference errors."""
    pass


# Module-level constants
SPECIES_NAME_MAPPING = {
    'human': 'homo-sapiens',
    'mouse': 'mus-musculus',
    'rat': 'rattus-norvegicus',
    'zebrafish': 'danio-rerio',
    'fly': 'drosophila-melanogaster',
    'worm': 'caenorhabditis-elegans',
    'yeast': 'saccharomyces-cerevisiae',
    'arabidopsis': 'arabidopsis-thaliana',
    'chicken': 'gallus-gallus',
    'cow': 'bos-taurus',
    'pig': 'sus-scrofa'
}

TIE_BREAKER_PRIORITY = ['yeast', 'worm', 'arabidopsis', 'fly', 'zebrafish', 'rat', 'mouse', 'chicken', 'cow', 'pig', 'human']

GENE_FORMAT_MAPPING = {
    'symbol': 'gene symbol',
    'ensembl': 'Ensembl Id',
    'entrez': 'Entrez Id'
}


def _select_best_species(species_scores: Dict[str, float]) -> str:
    """
    Select the best species from species scores using tie-breaking logic.
    
    Args:
        species_scores: Dictionary mapping species names to confidence scores
        
    Returns:
        Best species name (internal format, e.g., 'human')
    """
    if not species_scores:
        raise SpeciesInferenceError("No species scores provided")
    
    max_score = max(species_scores.values())
    best_candidates = [species for species, score in species_scores.items() if score == max_score]
    
    # Tie-breaking priority: prioritize species with more distinctive patterns
    best_species = best_candidates[0]  # Default to first candidate
    for priority_species in TIE_BREAKER_PRIORITY:
        if priority_species in best_candidates:
            best_species = priority_species
            break
    
    return best_species


class GeneIdentifierAnalyzer:
    """Analyzes gene identifiers to determine their type and infer species."""
    
    # Common gene identifier patterns
    ENSEMBL_PATTERNS = {
        'human': re.compile(r'^ENSG\d{11}$'),
        'mouse': re.compile(r'^ENSMUSG\d{11}$'),
        'rat': re.compile(r'^ENSRNOG\d{11}$'),
        'zebrafish': re.compile(r'^ENSDARG\d{11}$'),
        'fly': re.compile(r'^FBgn\d{7}$'),
        'worm': re.compile(r'^WBGene\d{8}$'),
        'yeast': re.compile(r'^Y[A-P][LR]\d{3}[CW]?$'),
        'arabidopsis': re.compile(r'^AT[1-5]G\d{5}$'),
        'chicken': re.compile(r'^ENSGALG\d{11}$'),
        'cow': re.compile(r'^ENSBTAG\d{11}$'),
        'pig': re.compile(r'^ENSSSCG\d{11}$')
    }
    
    # Entrez ID ranges for different species (approximate)
    ENTREZ_RANGES = {
        'human': (1, 100000),      # Human Entrez IDs typically 1-100k
        'mouse': (100000, 200000), # Mouse Entrez IDs typically 100k-200k  
        'rat': (200000, 300000),   # Rat Entrez IDs typically 200k-300k
        'zebrafish': (300000, 400000), # Zebrafish Entrez IDs typically 300k-400k
        'fly': (400000, 500000),   # Fly Entrez IDs typically 400k-500k
        'worm': (500000, 600000),  # Worm Entrez IDs typically 500k-600k
        'arabidopsis': (600000, 700000), # Arabidopsis Entrez IDs typically 600k-700k
        'chicken': (700000, 800000), # Chicken Entrez IDs typically 700k-800k
        'cow': (800000, 900000),   # Cow Entrez IDs typically 800k-900k
        'pig': (900000, 1000000),  # Pig Entrez IDs typically 900k-1M
        'yeast': (1000000, 1100000) # Yeast Entrez IDs typically 1M-1.1M
    }
    
    # Species-specific gene symbol patterns and known genes
    # Priority order: Most distinctive patterns first
    SPECIES_GENES = {
        'human': {
            'symbols': {
                'A3GALT2', 'ABCF2-H2BK1', 'ABHD14A-ACY1', 'ACOT1', 'ACSM2A', 'ACSM6', 'ACTR3C', 'ACY3', 'ADGRE1', 'ADGRE2',
                'ADGRG4', 'ADH1A', 'ADH1B', 'AGAP4', 'AGAP5', 'AGAP6', 'AGAP9', 'AK4P3', 'AKR1B15', 'AKR1C1',
                'AKR7A3', 'AKR7L', 'ALG10B', 'ALG1L2', 'ALPG', 'ALPP', 'AMY1B', 'AMY1C', 'ANKHD1-EIF4EBP3', 'ANKRD18A',
                'ANKRD18B', 'ANKRD20A1', 'ANKRD30A', 'ANKRD30B', 'ANKRD30BL', 'ANKRD36', 'ANKRD36B', 'ANKRD36C', 'ANKRD62', 'ANKRD7',
                'ANO7', 'ANP32D', 'ANXA2R', 'ANXA8', 'APOBEC3A', 'APOBEC3B', 'APOBEC3C', 'APOBEC3D', 'APOBEC3F', 'APOBEC3G',
                'APOBR', 'APOC1', 'APOC4-APOC2', 'APOL1', 'APOL2', 'APOL4', 'APOL5', 'APOL6', 'AQP12B', 'AQP7B',
                'ARGFX', 'ARHGAP11A-SCG5', 'ARHGAP19-SLIT1', 'ARHGEF35', 'ARL17A', 'ARL17B', 'ARL2-SNX15', 'ARL5C', 'ARMCX5-GPRASP2', 'ARMS2',
                'ARPC4-TTLL3', 'ARPIN-AP3S2', 'ARSF', 'ASAH2', 'ASAH2B', 'ASDURF', 'ATAD3C', 'ATF7-NPFF', 'ATP5MF-PTCD1', 'ATP5MGL',
                'ATP6V1G2-DDX39B', 'ATXN3L', 'BCL2L2-PABPN1', 'BIVM-ERCC5', 'BLID', 'BLOC1S5-TXNDC5', 'BOD1L2', 'BOLA2', 'BOLA2-SMG1P6', 'BORCS7-ASMT',
                'BORCS8-MEF2B', 'BPY2', 'BPY2B', 'BPY2C', 'BTN2A1', 'BTN3A1', 'BTNL3', 'BTNL8', 'BUB1B-PAK6', 'C10orf67',
                'C10orf95', 'C11orf21', 'C11orf96', 'C12orf42', 'C12orf71', 'C16orf54', 'C16orf95', 'C17orf50', 'C17orf99', 'C19orf18',
                'C19orf33', 'C19orf73', 'C1GALT1C1L', 'C1QTNF3-AMACR', 'C1orf162', 'C1orf202', 'C1orf226', 'C22orf42', 'C2CD4B', 'C2orf78',
                'C2orf92', 'C3orf20', 'C3orf22', 'C4B', 'C4B_2', 'C4orf50', 'C6orf118', 'C6orf141', 'C7orf33', 'C8orf33',
                'C8orf44-SGK3', 'C9orf153', 'C9orf43', 'C9orf50', 'C9orf57', 'CALCB', 'CALML6', 'CARD16', 'CARD18', 'CASP12',
                'CBLL2', 'CCDC144A', 'CCDC168', 'CCDC169-SOHLH2', 'CCDC17', 'CCDC187', 'CCDC200', 'CCDC201', 'CCDC7', 'CCDC74A',
                'CCDC74B', 'CCL13', 'CCL15', 'CCL15-CCL14', 'CCL18', 'CCL23', 'CCL3L3', 'CCL4L2', 'CCL7', 'CCNYL1B',
                'CCZ1B', 'CD1C', 'CD200R1L', 'CD300LD', 'CD300LF', 'CD33', 'CD8B2', 'CDC14C', 'CDK11A', 'CDK11B',
                'CDRT15', 'CDRT15L2', 'CDY1', 'CDY1B', 'CDY2A', 'CDY2B', 'CEACAM21', 'CEACAM3', 'CEACAM4', 'CEACAM5',
                'CEACAM6', 'CEACAM7', 'CEACAM8', 'CELA2B', 'CELA3A', 'CEMP1', 'CENPS-CORT', 'CENPVL1', 'CENPVL2', 'CENPVL3',
                'CEP295NL', 'CFAP144P1', 'CFAP298-TCP10L', 'CFAP47', 'CFAP97D2', 'CFC1B', 'CFHR1', 'CFHR2', 'CFHR3', 'CFHR4',
            },
            'patterns': [
                re.compile(r'^[A-Z]{2,6}\d*$'),  # Uppercase with numbers
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'mouse': {
            'symbols': {
                'Aadacl2fm1', 'Aadacl4fm4', 'Aadacl4fm5', 'Abca4', 'Abca8b', 'Abcb4', 'Abhd17a', 'Ablim1', 'Abo', 'Acad12',
                'Acat3', 'Acbd3', 'Accs', 'Acnat2', 'Acot10', 'Acot6', 'Actb', 'Acte1', 'Adam20', 'Adam26b',
                'Adam29', 'Adam3', 'Adam39', 'Adam4', 'Adam6a', 'Adam6b', 'Adck2', 'Adgrb1', 'Adgrg4', 'Adh1',
                'Adh6b', 'Adprhl1', 'Adrm1b', 'Afg1l', 'Agbl4', 'Ahcyl', 'Akr1c18', 'Akr1c20', 'Akr1c21', 'Akr1c6',
                'Akr1cl', 'Akr1e1', 'Akr7a5', 'Aldh3b3', 'Aldoart1', 'Alg10b', 'Allc', 'Alox8', 'Alppl2', 'Alyref2',
                'Alyreffm1', 'Alyreffm10', 'Alyreffm11', 'Alyreffm12', 'Alyreffm13', 'Alyreffm14', 'Alyreffm15', 'Alyreffm16', 'Alyreffm17', 'Alyreffm2',
                'Alyreffm3', 'Alyreffm4', 'Alyreffm5', 'Alyreffm6', 'Alyreffm7', 'Alyreffm8', 'Alyreffm9', 'Amd2', 'Amy2a1', 'Amy2a2',
                'Amy2a4', 'Amy2a5', 'Anapc15-ps', 'Ang4', 'Ang5', 'Ang6', 'Ankfn1', 'Ankrd48', 'Ankrd68', 'Ankrd69',
                'Anks1', 'Antxrl', 'Anxa2r1', 'Anxa2r2', 'Aoc1l1', 'Aoc1l2', 'Aoc1l3', 'Aoc2', 'Aoc3', 'Ap5z1',
                'Aph1c', 'Apoc2l', 'Apol10a', 'Apol10b', 'Apol11a', 'Apol11b', 'Apol6', 'Apol7a', 'Apol7b', 'Apol7c',
                'Apol7e', 'Apol8', 'Apol9b', 'Aqp12', 'Arf2', 'Arhgap33os', 'Arln', 'Art2a', 'Arxes1', 'Arxes2',
                'Ascl3', 'Asdurf', 'Asz1', 'Atg4a-ps', 'Atg5lrt', 'Atmin', 'Atp5mj', 'Atp6v0a2', 'Atp6v0e', 'Atp6v1f',
                'Aurkaip1', 'Bacc1', 'Bcl2a1a', 'Bcl2a1b', 'Bcl2a1c', 'Bcl2a1d', 'Bclaf1', 'Bco2', 'Bend2', 'Bex1',
                'Bex6', 'Bglap2', 'Bglap3', 'Bhlhb9', 'Bhmt1b', 'Bloc1s1', 'Bltp2', 'Bmp8a', 'Bmp8b', 'Bod1l',
                'Bpifb9b', 'Brd8dc', 'Btbd35f1', 'Btbd35f10', 'Btbd35f11', 'Btbd35f12', 'Btbd35f13', 'Btbd35f14', 'Btbd35f15', 'Btbd35f16',
                'Btbd35f17', 'Btbd35f18', 'Btbd35f19', 'Btbd35f2', 'Btbd35f20', 'Btbd35f21', 'Btbd35f22', 'Btbd35f23', 'Btbd35f24', 'Btbd35f25',
                'Btbd35f26', 'Btbd35f27', 'Btbd35f29', 'Btbd35f3', 'Btbd35f4', 'Btbd35f5', 'Btbd35f6', 'Btbd35f7', 'Btbd35f8', 'Btf3l4b',
                'Btnl1', 'Btnl6', 'C1ra', 'C1rb', 'C1s1', 'C1s2', 'C4a', 'C4bp', 'C5ar1', 'Calcoco2',
                'Calm4', 'Calm5', 'Caln1', 'Caml', 'Capns2', 'Capza1b', 'Casp16', 'Casp9', 'Catspere1', 'Catspere2',
                'Catsperg1', 'Catsperg2', 'Catsperz', 'Cbr1b', 'Cby1', 'Ccdc106', 'Ccdc121rt2', 'Ccdc121rt3', 'Ccdc15', 'Ccdc153',
            },
            'patterns': [
                re.compile(r'^[A-Z][a-z]+\d*$'),  # Mouse: Uppercase followed by lowercase
                re.compile(r'^[A-Z]{2,4}\d*[a-z]*$'),  # Mouse: Short mixed case
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'rat': {
            'symbols': {
                'A1bgl1', 'AAdacl4fm3l1', 'AKR1C3l1', 'Aadacl2fm2l1', 'Aadacl2fm3l1', 'Aadacl2l', 'Abca8', 'Abcg3l1', 'Abcg3l3', 'Abcg3l4',
                'Abo2', 'Abo3', 'Acat2l1', 'Actg1l1', 'Actl9b', 'Adam3a', 'Adam6', 'Adh1-ps1', 'Adh1c', 'Adh6',
                'Adnp2l1', 'Akr1b1-ps2', 'Akr1c1', 'Akr1c12l1', 'Akr1c15', 'Akr1c2', 'Akr1c3', 'Akr1e2', 'Akr7a2', 'Akr7a3',
                'Alg13l1', 'Alox15b', 'Alpg', 'Alpp', 'Andpro', 'Ankh', 'Ankra2l1', 'Ankrd36', 'Ankrd62l', 'Anks1a',
                'Anlnl1', 'Antxrl1', 'Antxrl3', 'Apbh', 'Apcdd1l', 'Apeg3', 'Apex2l1', 'Aph1bl2', 'Apol2', 'Apol3l1',
                'Apol7al1', 'Apol7bl1', 'Apool2', 'Aqp10', 'Aqp12a', 'Arhgap20l1', 'Arpc5l-ps1', 'Arsl', 'Asip', 'Asmtl',
                'Atp5mc2l1', 'Atp5mgl3', 'Atp5mj-ps1', 'Atp5mk-ps2', 'Atp5mpl', 'Atp6v0e1', 'Bcl2a1', 'Bin2a', 'Bles03', 'Bmt2',
                'Bod1l1', 'Bola2-ps1', 'Bpifa2f', 'Btf3-ps1', 'Btf3l6', 'Btg1b', 'Btg1c', 'Btnl3', 'Btnl5', 'Btnl8',
                'Bves', 'C10h1orf35', 'C10h5orf15', 'C10h5orf47', 'C10h5orf58', 'C11h3orf38', 'C11h3orf52', 'C11h3orf70', 'C12h7orf50', 'C13h1orf21',
                'C13h1orf53', 'C13h1orf74', 'C13h2orf76', 'C14h2orf73', 'C14h2orf74', 'C14h4orf19', 'C14h4orf36', 'C14h4orf50', 'C14h5orf52', 'C14h7orf57',
                'C15h3orf49', 'C15h8orf58', 'C15h8orf74', 'C16h8orf48', 'C17h5orf24', 'C17h6orf52', 'C17h6orf62', 'C17h7orf25', 'C18h5orf46', 'C18h5orf63',
                'C19h4orf51', 'C1h10orf88', 'C1h10orf90', 'C1h10orf95', 'C1h11orf16', 'C1h11orf24', 'C1h11orf42', 'C1h11orf58', 'C1h11orf86', 'C1h11orf98',
                'C1h15orf40-ps2', 'C1h16orf54', 'C1h16orf82', 'C1h19orf12', 'C1h19orf18', 'C1h19orf33', 'C1h19orf47', 'C1h19orf81', 'C1h19orf84', 'C1h19orf85',
                'C1h2orf78', 'C1h2orf78l1', 'C1h6orf120', 'C1h9orf40', 'C1h9orf85', 'C1r', 'C1s', 'C20h6orf89', 'C2h1orf162', 'C2h1orf52',
                'C2h1orf54', 'C2h1orf56', 'C2h3orf33', 'C2h3orf80', 'C2h4orf17', 'C2h4orf3', 'C2h4orf33', 'C2h4orf46', 'C2h4orf54', 'C2h5orf22',
                'C2h5orf34', 'C3h11orf91', 'C3h11orf96', 'C3h15orf48', 'C3h15orf62', 'C3h20orf96', 'C3h9orf50', 'C3h9orf78', 'C4bpa', 'C4bpb',
                'C4h12orf60', 'C4h12orf71', 'C4h2orf42', 'C4h2orf68', 'C4h2orf81', 'C4h3orf20', 'C4h3orf22', 'C5h1orf127', 'C5h1orf141', 'C5h1orf159',
                'C5h1orf174', 'C5h1orf185', 'C5h1orf210', 'C5h1orf216', 'C5h1orf232', 'C5h1orf50', 'C5h1orf87', 'C5h1orf94', 'C5h6orf163', 'C5h8orf34',
                'C5h8orf88', 'C5h8orf89', 'C5h9orf152', 'C5h9orf43', 'C5l1', 'C6h14orf28', 'C6h2orf50', 'C7h12orf42', 'C7h12orf50', 'C7h12orf56',
                'C7h12orf75', 'C7h19orf25', 'C7h22orf23', 'C7h8orf33', 'C7h8orf76', 'C7h8orf82', 'C8h11orf54', 'C8h11orf65', 'C8h11orf71', 'C8h11orf87',
            },
            'patterns': [
                re.compile(r'^[A-Z][a-z]+\d*$'),  # Rat: Similar to mouse
                re.compile(r'^[A-Z]{2,4}\d*[a-z]*$'),  # Rat: Short mixed case
                re.compile(r'^RGD\d*$'),  # Rat-specific RGD pattern
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'zebrafish': {
            'symbols': {
                'C10H5orf63', 'C11H1orf50', 'C13H9orf72', 'C14HXorf65', 'C16H1orf43', 'C18H3orf33', 'C19H6orf62', 'C1H4orf33', 'C20H6orf58', 'C21H5orf24',
                'C22H1orf53', 'C23H1orf52', 'C24H8orf34', 'C2H1orf35', 'C2H5orf22', 'C2H7orf25', 'C3H16orf89', 'C3H17orf75', 'C5H11orf54', 'C5H2orf42',
                'C5H2orf68', 'C5H9orf78', 'C5H9orf85', 'C6H3orf18', 'C7H11orf68', 'C7H16orf87', 'C7H4orf54', 'C9H2orf49', 'C9H2orf69', 'C9H2orf76',
                'C9H3orf38', 'C9H3orf70', 'a1cf', 'a2ml', 'a2ml2', 'aaas', 'aacs', 'aadac', 'aadacl4', 'aadat',
                'aagab', 'aak1a', 'aak1b', 'aamdc', 'aamp', 'aanat1', 'aanat2', 'aar2', 'aars1', 'aars2',
                'aarsd1', 'aasdh', 'aasdhppt', 'aass', 'aatf', 'aatka', 'aatkb', 'abat', 'abca12', 'abca1a',
                'abca1b', 'abca2', 'abca3b', 'abca4a', 'abca4b', 'abca5', 'abca7', 'abcb10', 'abcb11a', 'abcb11b',
                'abcb4', 'abcb5', 'abcb6a', 'abcb6b', 'abcb7', 'abcb8', 'abcb9', 'abcc1', 'abcc10', 'abcc12',
                'abcc13', 'abcc2', 'abcc3', 'abcc4', 'abcc5', 'abcc6a', 'abcc6b.1', 'abcc6b.2', 'abcc8', 'abcc8b',
                'abcc9', 'abcd1', 'abcd2', 'abcd3a', 'abcd3b', 'abcd4', 'abce1', 'abcf1', 'abcf2a', 'abcf2b',
                'abcf3', 'abcg1', 'abcg2a', 'abcg2b', 'abcg2c', 'abcg2d', 'abcg4a', 'abcg4b', 'abcg5', 'abcg8',
                'abch1', 'abhd10a', 'abhd11', 'abhd12', 'abhd13', 'abhd14a', 'abhd14b', 'abhd15a', 'abhd16a', 'abhd17aa',
                'abhd17ab', 'abhd17b', 'abhd17c', 'abhd18', 'abhd2a', 'abhd2b', 'abhd3', 'abhd4', 'abhd5b', 'abhd6a',
                'abhd6b', 'abhd8a', 'abhd8b', 'abi1a', 'abi1b', 'abi2a', 'abi2b', 'abi3a', 'abi3b', 'abi3bpa',
                'abi3bpb', 'abitram', 'abl1', 'abl2', 'ablim1a', 'ablim1b', 'ablim2', 'ablim3', 'abr', 'abraa',
                'abrab', 'abracl', 'abraxas1', 'abraxas2', 'abt1', 'abtb1', 'abtb2a', 'abtb2b', 'acaa1', 'acaa2',
                'acaca', 'acacb', 'acad11', 'acad8', 'acad9', 'acadl', 'acadm', 'acads', 'acadsb', 'acadvl',
                'acana', 'acanb', 'acap1', 'acap2a', 'acap2b', 'acap3a', 'acap3b', 'acat1', 'acat2', 'acbd3',
                'acbd4', 'acbd5a', 'acbd5b', 'acbd6', 'acbd7', 'accs', 'acd', 'ace', 'ace2', 'acer1',
                'acer2', 'acer3', 'ache', 'acin1a', 'acin1b', 'ackr3a', 'ackr3b', 'ackr4a', 'ackr4b', 'aclya',
            },
            'patterns': [
                re.compile(r'^[A-Z][a-z]+\d*[a-z]*$'),  # Zebrafish: Mixed case with numbers
                re.compile(r'^[A-Z]{2,4}\d*[a-z]*$'),  # Zebrafish: Short uppercase with numbers
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'fly': {
            'symbols': {
                '14-3-3epsilon', '14-3-3zeta', '2mit', '4E-T', '5-HT1A', '5-HT1B', '5-HT2A', '5-HT2B', '5-HT7', '5PtaseI',
                '825-Oak', 'AANATL2', 'AANATL3', 'AANATL4', 'AANATL5', 'AANATL6', 'AANATL7', 'ACC', 'ACXA', 'ACXB',
                'ACXC', 'ACXD', 'ACXE', 'AGBE', 'AIF', 'AIMP3', 'ALiX', 'AMPKalpha', 'AMPdeam', 'AOX3',
                'AP-1-2beta', 'AP-1gamma', 'AP-1mu', 'AP-1sigma', 'AP-2alpha', 'AP-2mu', 'AP-2sigma', 'APP-BP1', 'AQP', 'ARY',
                'ASPP', 'ATP8A', 'ATP8B', 'ATPCL', 'ATPsynB', 'ATPsynC', 'ATPsynCF6', 'ATPsynCF6L', 'ATPsynD', 'ATPsynE',
                'ATPsynF', 'ATPsynG', 'ATPsynGL', 'ATPsynO', 'ATPsynbeta', 'ATPsynbetaL', 'ATPsyndelta', 'ATPsynepsilonL', 'ATPsyngamma', 'ATbp',
                'Abd-B', 'Abi', 'Abl', 'Abp1', 'Ac13E', 'Ac3', 'Ac76E', 'Ac78C', 'AcCoAS', 'Acam',
                'Acbp1', 'Acbp2', 'Acbp3', 'Acbp4', 'Acbp5', 'Acbp6', 'Acer', 'Acf', 'Achl', 'Ack',
                'Ack-like', 'Acn', 'Acox57D-d', 'Acox57D-p', 'Acp24A4', 'Acp26Aa', 'Acp26Ab', 'Acp29AB', 'Acp32CD', 'Acp33A',
                'Acp36DE', 'Acp53C14a', 'Acp53C14b', 'Acp53C14c', 'Acp53Ea', 'Acp54A1', 'Acp62F', 'Acp63F', 'Acp65Aa', 'Acp76A',
                'Acp95EF', 'Acp98AB', 'Acph-1', 'Acsl', 'Acsx1L', 'Acsx1R', 'Acsx2', 'Acsx3', 'Acsx4', 'Act42A',
                'Act57B', 'Act5C', 'Act79B', 'Act87E', 'Act88F', 'Actbeta', 'Actn', 'Acyp', 'Ada1-1', 'Ada1-2',
                'Ada2a', 'Ada2b', 'Ada3', 'AdamTS-A', 'AdamTS-B', 'Adf1', 'Adgf-A', 'Adgf-A2', 'Adgf-B', 'Adgf-C',
                'Adgf-D', 'Adgf-E', 'Adh', 'Adhr', 'AdipoR', 'Adk1', 'Adk2', 'Adk3', 'AdoR', 'Adss',
                'Aduk', 'Ae2', 'Aef1', 'Afti', 'Ag5r', 'Ag5r2', 'AgmNAT', 'AhcyL1', 'AhcyL2', 'Akap200',
                'Akh', 'AkhR', 'Akr1B', 'Akt', 'AlaRS', 'AlaRS-m', 'Aladin', 'Alas', 'Ald1', 'Ald2',
                'Aldh', 'Aldh-III', 'Aldh7A1', 'Alg-2', 'Alg7', 'Alh', 'AlkB', 'Alms1a', 'Alms1b', 'Alp1',
                'Alp10', 'Alp11', 'Alp12', 'Alp13', 'Alp2', 'Alp4', 'Alp5', 'Alp6', 'Alp7', 'Alp8',
                'Alp9', 'Alr', 'Ama', 'Amnionless', 'Amun', 'Amy-d', 'Amy-p', 'Amyrel', 'Ance', 'Ance-2',
                'Ance-3', 'Ance-4', 'Ance-5', 'Andorra', 'Anp', 'Ant2', 'Antdh', 'Antp', 'AnxB10', 'AnxB11',
            },
            'patterns': [
                re.compile(r'^[a-z]+(-[a-z]+)*$'),  # Fly pattern: lowercase with optional hyphens
                re.compile(r'^[A-Z][a-z]+$'),  # Fly pattern: single uppercase followed by lowercase
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'worm': {
            'symbols': {
                'AC3.12', 'AC3.5', 'AC7.3', 'AC8.10', 'AC8.3', 'AC8.4', 'AH10.2', 'AH10.4', 'AH6.17', 'AH6.3',
                'AH9.1', 'AH9.3', 'AH9.4', 'AH9.6', 'BE10.1', 'BE10.3', 'BE10.4', 'BE10.5', 'C26H9A.3', 'C26H9A.4',
                'CC8.2', 'CC8.3', 'CD4.1', 'CD4.11', 'CD4.8', 'DC2.5', 'DH11.2', 'DH11.5', 'DY3.8', 'EEED8.12',
                'EEED8.13', 'EEED8.14', 'EEED8.15', 'EEED8.18', 'EEED8.2', 'EEED8.3', 'EEED8.4', 'EGAP1.1', 'EGAP1.5', 'EGAP2.1',
                'EGAP2.2', 'EGAP4.1', 'EGAP5.1', 'EGAP9.3', 'EGAP9.4', 'JC8.16', 'JC8.4', 'JC8.7', 'M6.4', 'M7.7',
                'M7.8', 'PDB1.1', 'R05G9R.1', 'VC5.2', 'VF13D12L.3', 'VH15N14R.1', 'VM106R.1', 'VW02B12L.2', 'VY10G11R.1', 'VZK822L.2',
                'Y108G3AL.2', 'Y110A2AL.1', 'Y110A2AL.2', 'Y110A2AL.3', 'Y110A2AL.4', 'Y110A2AL.5', 'Y110A2AL.6', 'Y110A2AL.7', 'Y110A2AL.9', 'Y110A2AR.1',
                'Y11D7A.3', 'Y11D7A.5', 'Y11D7A.7', 'Y11D7A.8', 'Y11D7A.9', 'Y12A6A.1', 'Y13C8A.1', 'Y15E3A.3', 'Y15E3A.4', 'Y15E3A.5',
                'Y16B4A.2', 'Y17D7B.2', 'Y17D7B.3', 'Y17D7B.4', 'Y17D7C.1', 'Y17D7C.2', 'Y17D7C.6', 'Y17G7B.3', 'Y17G7B.8', 'Y17G9A.2',
                'Y17G9A.3', 'Y17G9A.4', 'Y17G9B.2', 'Y17G9B.4', 'Y17G9B.5', 'Y17G9B.8', 'Y18H1A.2', 'Y18H1A.4', 'Y18H1A.8', 'Y18H1A.9',
                'Y19D2B.1', 'Y19D2B.2', 'Y1A5A.1', 'Y1B5A.1', 'Y20C6A.1', 'Y20C6A.4', 'Y22D7AL.1', 'Y22D7AL.10', 'Y22D7AL.11', 'Y22D7AL.15',
                'Y22D7AL.16', 'Y22D7AL.3', 'Y22D7AL.4', 'Y22D7AL.6', 'Y22D7AL.7', 'Y22D7AL.9', 'Y22D7AR.14', 'Y22D7AR.2', 'Y22D7AR.3', 'Y22D7AR.6',
                'Y22D7AR.7', 'Y23B4A.1', 'Y23H5A.2', 'Y23H5A.8', 'Y23H5B.1', 'Y23H5B.3', 'Y23H5B.5', 'Y23H5B.8', 'Y24D9A.6', 'Y24D9A.7',
                'Y24D9B.1', 'Y25C1A.2', 'Y25C1A.6', 'Y25C1A.7', 'Y25C1A.8', 'Y26D4A.3', 'Y26D4A.5', 'Y26E6A.2', 'Y26E6A.3', 'Y27F2A.6',
                'Y27F2A.8', 'Y27F2A.9', 'Y2H9A.4', 'Y2H9A.6', 'Y32F6A.4', 'Y32F6A.5', 'Y32F6B.1', 'Y32F6B.9', 'Y32G9A.2', 'Y32G9A.3',
                'Y32G9A.5', 'Y32G9B.1', 'Y34B4A.2', 'Y34B4A.4', 'Y34B4A.7', 'Y34D9A.3', 'Y34D9A.7', 'Y34D9A.8', 'Y36E3A.2', 'Y37A1A.2',
                'Y37A1A.3', 'Y37A1A.4', 'Y37A1B.4', 'Y37A1B.7', 'Y37D8A.2', 'Y37D8A.3', 'Y37D8A.4', 'Y37D8A.5', 'Y37D8A.6', 'Y37D8A.8',
                'Y37E11AL.1', 'Y37E11AL.2', 'Y37E11AL.3', 'Y37E11AL.9', 'Y37E11AM.2', 'Y37E11AM.3', 'Y37E11AM.4', 'Y37E11AR.7', 'Y37H2A.1', 'Y37H2A.7',
                'Y37H2C.4', 'Y37H9A.1', 'Y37H9A.5', 'Y38C1AA.1', 'Y38C1AA.12', 'Y38C1AA.19', 'Y38C1AA.6', 'Y38C1AA.7', 'Y38C1AB.3', 'Y38C1AB.5',
                'Y38C1BA.1', 'Y38C1BA.4', 'Y38F1A.1', 'Y38F1A.2', 'Y38F1A.4', 'Y38F1A.7', 'Y38F1A.8', 'Y38F2AL.12', 'Y38F2AR.10', 'Y38F2AR.12',
            },
            'patterns': [
                re.compile(r'^[a-z]+-\d+$'),  # Worm pattern: lowercase-hyphen-number
                re.compile(r'^[a-z]+-\d+\.\d+$'),  # Worm pattern with decimals
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'yeast': {
            'symbols': {
                'AAD10', 'AAD14', 'AAD15', 'AAD3', 'AAD4', 'AAH1', 'AAT2', 'ABD1', 'ABM1', 'ABP140',
                'ABZ1', 'ABZ2', 'ACB1', 'ACF2', 'ACF4', 'ACH1', 'ACK1', 'ACL4', 'ACM1', 'ADD37',
                'ADD66', 'ADE1', 'ADE12', 'ADE13', 'ADE16', 'ADE17', 'ADE2', 'ADE3', 'ADE4', 'ADE5,7',
                'ADE6', 'ADE8', 'ADH2', 'ADH3', 'ADO1', 'ADP1', 'ADY2', 'ADY3', 'ADY4', 'AEP1',
                'AEP2', 'AEP3', 'AFB1', 'AFG1', 'AFG2', 'AFG3', 'AFI1', 'AFR1', 'AFT1', 'AFT2',
                'AGA1', 'AGA2', 'AGC1', 'AGE1', 'AGE2', 'AGX1', 'AHA1', 'AHC1', 'AHC2', 'AHK1',
                'AHT1', 'AI1', 'AI2', 'AI3', 'AI4', 'AI5_ALPHA', 'AI5_BETA', 'AIM10', 'AIM11', 'AIM14',
                'AIM17', 'AIM18', 'AIM19', 'AIM20', 'AIM21', 'AIM22', 'AIM23', 'AIM24', 'AIM25', 'AIM26',
                'AIM29', 'AIM3', 'AIM32', 'AIM33', 'AIM34', 'AIM36', 'AIM39', 'AIM4', 'AIM41', 'AIM44',
                'AIM45', 'AIM46', 'AIM6', 'AIM7', 'AIM9', 'AIP1', 'AIR2', 'AKL1', 'AKR1', 'ALD2',
                'ALD3', 'ALD4', 'ALD5', 'ALD6', 'ALF1', 'ALG7', 'ALK1', 'ALK2', 'ALO1', 'ALP1',
                'ALR1', 'ALR2', 'ALT1', 'ALT2', 'AMA1', 'AMD2', 'AME1', 'AMF1', 'AMS1', 'ANB1',
                'ANP1', 'ANR2', 'ANS1', 'ANY1', 'AOS1', 'APA2', 'APC5', 'APC9', 'APE3', 'APE4',
                'API2', 'APJ1', 'APL5', 'APL6', 'APM2', 'APM3', 'APM4', 'APN1', 'APN2', 'APQ12',
                'APQ13', 'AQR1', 'AQY1', 'AQY2', 'AQY3', 'ARA2', 'ARB1', 'ARC1', 'ARC15', 'ARC18',
                'ARC19', 'ARC35', 'ARC40', 'ARE1', 'ARE2', 'ARG3', 'ARG4', 'ARG5,6', 'ARG7', 'ARG8',
                'ARG80', 'ARG81', 'ARG82', 'ARH1', 'ARK1', 'ARN1', 'ARN2', 'ARO10', 'ARO7', 'ARO8',
                'ARO80', 'ARO9', 'ARP10', 'ARP3', 'ARR1', 'ARR2', 'ART10', 'ARX1', 'ASC1', 'ASF1',
                'ASF2', 'ASH1', 'ASI1', 'ASI2', 'ASI3', 'ASK1', 'ASK10', 'ASM4', 'ASP3-1', 'ASP3-2',
                'ASP3-3', 'ASP3-4', 'ASR1', 'AST1', 'AST2', 'ATC1', 'ATG1', 'ATG11', 'ATG15', 'ATG16',
            },
            'patterns': [
                re.compile(r'^Y[A-P][LR]\d{3}[CW]?(-[A-Z])?$'),  # Yeast systematic names
                re.compile(r'^[A-Z]{3,4}\d+$'),  # Yeast gene symbols
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'arabidopsis': {
            'symbols': {
                '2-Cys Prx B', '3BETAHSD/D1', '3BETAHSD/D2', '3xHMG-box1', '3xHMG-box2', '4CL1', '4CL2', '4CL3', '4CL5', '4CL8',
                '5-FCL', '5PTASE11', '5PTASE13', '5PTase12', '5PTase14', 'A/N-InvA', 'A/N-InvB', 'A/N-InvC', 'A/N-InvD', 'A/N-InvF',
                'AAC2', 'AACT1', 'AAE1', 'AAE11', 'AAE12', 'AAE13', 'AAE14', 'AAE15', 'AAE16', 'AAE17',
                'AAE18', 'AAE2', 'AAE3', 'AAE5', 'AAE7', 'AAH', 'AAK6', 'AAO2', 'AAO3', 'AAP2',
                'AAP3', 'AAP4', 'AAP5', 'AAP6', 'AAP7', 'AAP8', 'AAPT1', 'AAPT2', 'AAR3', 'AARE',
                'AAS', 'AAT', 'AATP1', 'ABA1', 'ABA2', 'ABA3', 'ABA4', 'ABAP1', 'ABC1', 'ABC4',
                'ABCA11', 'ABCB12', 'ABCB13', 'ABCB14', 'ABCB15', 'ABCB16', 'ABCB17', 'ABCB18', 'ABCB19', 'ABCB2',
                'ABCB20', 'ABCB21', 'ABCB22', 'ABCB23', 'ABCB24', 'ABCB25', 'ABCB26', 'ABCB27', 'ABCB28', 'ABCB29',
                'ABCB3', 'ABCC14', 'ABCC15', 'ABCC7', 'ABCE2', 'ABCE3', 'ABCF4', 'ABCF5', 'ABCG10', 'ABCG11',
                'ABCG12', 'ABCG13', 'ABCG14', 'ABCG15', 'ABCG16', 'ABCG17', 'ABCG18', 'ABCG19', 'ABCG20', 'ABCG21',
                'ABCG22', 'ABCG23', 'ABCG24', 'ABCG25', 'ABCG26', 'ABCG28', 'ABCG29', 'ABCG3', 'ABCG30', 'ABCG31',
                'ABCG32', 'ABCG33', 'ABCG34', 'ABCG35', 'ABCG37', 'ABCG38', 'ABCG39', 'ABCG40', 'ABCG41', 'ABCG42',
                'ABCG43', 'ABCG6', 'ABCG7', 'ABCG9', 'ABCI1', 'ABCI10', 'ABCI11', 'ABCI12', 'ABCI13', 'ABCI14',
                'ABCI15', 'ABCI16', 'ABCI17', 'ABCI18', 'ABCI19', 'ABCI20', 'ABCI21', 'ABCI4', 'ABCI5', 'ABCI7',
                'ABCI8', 'ABCI9', 'ABF3', 'ABF4', 'ABH1', 'ABI4', 'ABI5', 'ABIL1', 'ABIL2', 'ABIL3',
                'ABIL4', 'ABO1', 'ABO3', 'ABO5', 'ABO6', 'ABR1', 'ABS2', 'ACA.l', 'ACA10', 'ACA11',
                'ACA2', 'ACA3', 'ACA4', 'ACA5', 'ACA6', 'ACA7', 'ACA8', 'ACA9', 'ACBP1', 'ACBP2',
                'ACBP3', 'ACBP4', 'ACBP5', 'ACBP6', 'ACC2', 'ACD1', 'ACD1-LIKE', 'ACD11', 'ACD2', 'ACD32.1',
                'ACD5', 'ACD6', 'ACDO1', 'ACHT1', 'ACHT2', 'ACHT3', 'ACHT4', 'ACHT5', 'ACI1', 'ACL',
                'ACL5', 'ACLA-1', 'ACLA-2', 'ACLA-3', 'ACLB-1', 'ACLB-2', 'ACO3', 'ACOS5', 'ACR1', 'ACR10',
            },
            'patterns': [
                re.compile(r'^AT[1-5]G\d{5}$'),  # Arabidopsis systematic names
                re.compile(r'^[A-Z]{2,4}\d*$'),  # Arabidopsis gene symbols
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'chicken': {
            'symbols': {
                'A2ML2', 'A2ML3', 'A2ML4', 'AADACL3B', 'AADACL3C', 'AADACL4A', 'AADACL4B', 'AAED1', 'ABCB1LA', 'ACOD1L',
                'ACPP', 'ACRL', 'ACSS1A', 'ACSS1B', 'ACTG1L', 'ACTR10L', 'ADAM20L', 'ADAM32L1', 'ADAM32L2', 'ADAM9L',
                'ADAMTSL2L', 'ADAP1L', 'ADMP', 'ADPRHL', 'ADPRHL2', 'AES', 'AHR1A', 'AHR1B', 'AHR2', 'AIM1L',
                'ALPPL', 'ALS2L', 'ANGPT1L', 'ANKFN1L', 'ANKRD9L', 'APH1AL', 'APLNR2', 'APOA1BP', 'APOPT1', 'APOV1',
                'ARHGAP20B', 'ARHGAP33L1', 'ARHGAP33L4', 'ARHGAP33L6', 'ARHGAP39L', 'ARHGEF7L', 'ARL2BPL', 'ARL8BL', 'ARMC4', 'ART1A',
                'ART1B', 'ART1D', 'ART1E', 'ART1F', 'ATP5B', 'ATP5C1', 'ATP5D', 'ATP5E', 'ATP5F1AW', 'ATP5F1AZ',
                'ATP5G1', 'ATP5G2', 'ATP5G3', 'ATP5H', 'ATP5I', 'ATP5J', 'ATP5J2', 'ATP5L', 'ATP5O', 'ATP5S',
                'ATP5SL', 'ATPIF1', 'ATVR1', 'AVD', 'AVDL', 'AVPR2C', 'AVR2', 'AVR7', 'AvBD1', 'AvBD10',
                'AvBD11', 'AvBD13', 'AvBD14', 'AvBD3', 'AvBD4', 'AvBD6', 'AvBD7', 'AvBD8', 'AvBD9', 'B3GAT1L',
                'B3GNT10', 'B4GALTL', 'BDH1A', 'BDH1B', 'BF1', 'BF2', 'BG8', 'BGV', 'BIRC8', 'BKJ',
                'BLB1', 'BLB2', 'BLEC1', 'BLEC2', 'BLEC3', 'BPIFCB', 'BPIL3', 'BRE', 'BTBD11', 'BTG1L',
                'BTN1', 'BlSK1', 'C11orf1', 'C11orf49', 'C11orf57', 'C11orf63', 'C11orf70', 'C11orf74', 'C11orf88', 'C11orf94',
                'C12orf10', 'C12orf29', 'C12orf4', 'C12orf45', 'C12orf49', 'C12orf65', 'C12orf73', 'C14ORF166B', 'C14ORF39', 'C14orf2',
                'C15ORF40', 'C15orf52', 'C16orf52', 'C16orf62', 'C16orf72', 'C17ORF80', 'C18orf8', 'C19orf24', 'C19orf60', 'C19orf70',
                'C19orf71', 'C1H3ORF52', 'C1HXORF59', 'C1orf101', 'C1orf106', 'C1orf112', 'C1orf131', 'C1orf158', 'C1orf228', 'C20orf85',
                'C21orf59', 'C21orf62', 'C2H5ORF22', 'C2H6orf52', 'C2H8ORF46', 'C2H8ORF76', 'C2H8orf82', 'C2H9ORF30', 'C2orf50', 'C2orf70',
                'C2orf71', 'C3AR1L', 'C3H1ORF115', 'C3H8ORF80', 'C3orf14', 'C3orf67', 'C4BPG', 'C4BPM', 'C4BPS', 'C4H2orf81',
                'C4H4ORF46', 'C4H4ORF50', 'C4H4orf54', 'C4ORF19', 'C4orf45', 'C4orf47', 'C4orf48', 'C5H11ORF9', 'C5H14orf4', 'C5orf49',
                'C5orf51', 'C6orf106', 'C6orf203', 'C7orf26', 'C7orf31', 'C7orf50', 'C7orf62', 'C7orf72', 'C7orf73', 'C8orf22',
                'C8orf37', 'C8orf4', 'C8orf59', 'C9ORF152', 'C9ORF58', 'C9orf116', 'C9orf142', 'C9orf24', 'C9orf64', 'C9orf69',
            },
            'patterns': [
                re.compile(r'^[A-Z]{2,6}\d*$'),  # Uppercase with numbers
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'cow': {
            'symbols': {
                '20ALPHA-HSD', 'ABCA14', 'ABCA16', 'ACTE1', 'ADAM1A', 'ADGRF2', 'ADID', 'ALOX12E', 'ANG2', 'APOBEC3Z1',
                'AQP12', 'BCAR4', 'BCNT2', 'BDA20', 'BOLA-DMB', 'BOLA-DOA', 'BOLA-DOB', 'BOLA-DQA1', 'BOLA-DQA5', 'BOLA-DQB',
                'BOLA-DQB1', 'BOLA-DRA', 'BOLA-DYA', 'BOLA-DYB', 'BOSTAUV1R416', 'BOTA-T2R10B', 'BPIFA2A', 'BPIFA2B', 'BPIFA2C', 'BPIFB5',
                'BRB', 'BREH1', 'BSP3', 'BSP5', 'BTY3', 'BTY4', 'BTY5', 'BTY8', 'BTY9', 'BoLA',
                'C11H2orf15', 'C11H2orf42', 'C11H2orf49', 'C11H2orf50', 'C11H2orf68', 'C11H2orf73', 'C11H2orf74', 'C11H2orf78', 'C11H2orf81', 'C11H9orf50',
                'C11H9orf78', 'C14H8orf33', 'C14H8orf34', 'C14H8orf76', 'C14H8orf82', 'C14H8orf88', 'C14H8orf89', 'C16H1orf21', 'C16H1orf53', 'C16H1orf74',
                'C17H4orf33', 'C17H4orf46', 'C17H5orf52', 'C1H21orf91', 'C1H3orf33', 'C1H3orf38', 'C1H3orf52', 'C1H3orf70', 'C1H3orf80', 'C20H5orf22',
                'C20H5orf34', 'C20H5orf47', 'C22H3orf18', 'C22H3orf20', 'C22H3orf22', 'C22H3orf49', 'C22H3orf62', 'C22H3orf84', 'C23H6orf15', 'C23H6orf47',
                'C23H6orf62', 'C23H6orf89', 'C25H7orf50', 'C27H8orf48', 'C2H2orf69', 'C2H2orf72', 'C2H2orf76', 'C2H2orf80', 'C2H2orf88', 'C3H1orf122',
                'C3H1orf141', 'C3H1orf146', 'C3H1orf185', 'C3H1orf210', 'C3H1orf216', 'C3H1orf43', 'C3H1orf50', 'C3H1orf52', 'C3H1orf54', 'C3H1orf56',
                'C3H1orf87', 'C3H1orf94', 'C4H7orf25', 'C4H7orf57', 'C5H12orf50', 'C5H12orf54', 'C5H12orf56', 'C5H12orf57', 'C5H12orf60', 'C5H12orf71',
                'C5H12orf75', 'C5H22orf23', 'C6H4orf17', 'C6H4orf19', 'C6H4orf3', 'C7H19orf25', 'C7H19orf38', 'C7H19orf44', 'C7H19orf67', 'C7H1orf35',
                'C7H5orf15', 'C7H5orf24', 'C7H5orf46', 'C7H5orf63', 'C8H8orf58', 'C8H8orf74', 'C8H9orf152', 'C8H9orf40', 'C8H9orf43', 'C8H9orf85',
                'C9H6orf118', 'C9H6orf120', 'C9H6orf163', 'C9H6orf58', 'CALM', 'CASP16', 'CATHL1', 'CATHL2', 'CATHL3', 'CATHL4',
                'CATHL5', 'CATHL6', 'CCDC105', 'CCDC115', 'CD1B5', 'CDK11', 'CFDP2', 'CGN1', 'CL43', 'CL46',
                'CLCA3', 'CLECL1', 'COX8B', 'CRYGF', 'CSN1S2', 'CSNK1B', 'CUPIN1', 'CXHXorf38', 'CXHXorf58', 'CXHXorf65',
                'CXHXorf66', 'CYCT', 'CYHR1', 'CYM', 'CYP21A1P', 'CYP2B39', 'CYP2C23', 'CYP2C282', 'CYP2C87', 'CYP2C88',
                'CYP2C89', 'CYP2C90', 'CYP2D43', 'CYP2J30', 'CYP3A28', 'CYP3A74', 'CYP4A59', 'DBIL5', 'DEFB', 'DEFB10',
                'DEFB122', 'DEFB122A', 'DEFB13', 'DEFB7', 'DOC2G', 'DSB', 'EBD', 'FADS2B', 'FAM205A', 'FAM205C',
                'GAT', 'GGNBP1', 'GHE4', 'GLYCAM1', 'GPB5', 'GRO1', 'GULO', 'H1-9', 'H2AC10', 'H2B',
            },
            'patterns': [
                re.compile(r'^[A-Z]{2,6}\d*$'),  # Uppercase with numbers
            ],
            'priority': 1  # High priority for species-specific genes
        },
        'pig': {
            'symbols': {
                'ABCG2A', 'AMCF-II', 'ARSE', 'AWN', 'C8H4orf50', 'CBR2', 'CCL3L1', 'COX8H', 'CPHL1', 'CYP2A19',
                'CYP2B6B', 'CYP2C42', 'CYP2C49', 'CYP2J34', 'CYP3A29', 'CYP4A21', 'CYP4A24', 'FOLH1B', 'GKN3', 'GP91-PHOX',
                'HSP70.2', 'IFN-ALPHA-1', 'IFN-ALPHA-10', 'IFN-ALPHA-13', 'IFN-ALPHA-15', 'IFN-ALPHA-17', 'IFN-ALPHA-4', 'IFN-ALPHA-8', 'IFN-ALPHA-9', 'IFN-ALPHAOMEGA',
                'IFN-DELTA-1', 'IFN-DELTA-10', 'IFN-DELTA-11', 'IFN-DELTA-4', 'IFN-DELTA-5', 'IFN-DELTA-6', 'IFN-DELTA-7', 'IFN-DELTA-8', 'IFN-DELTA-9', 'IFN-OMEGA-2',
                'IFN-OMEGA-3', 'IFN-OMEGA-5', 'IFN-OMEGA-6', 'IFN-OMEGA-7', 'IGHM', 'IL1B2', 'IL28B', 'IL29', 'ISG12(A)', 'MOBKL3',
                'NPG3-PG-2', 'OLF42-1', 'OLF42-3', 'PGB', 'PHEROC', 'PMAP-23', 'PMAP-36', 'PPAG3', 'PSP-I', 'PSP-II',
                'SERPINA3-2', 'SHAS2', 'SLA-DMA', 'SLA-DMB', 'SLA-DOA', 'SLA-DQB1', 'SPAI-2', 'SPMI', 'SPRP', 'TP23',
                'TRAV4', 'TRBV19', 'TRBV27', 'TRBV30', 'TRBV9', 'UABP-2', 'UGT2B31',
            },
            'patterns': [
                re.compile(r'^[A-Z]{2,6}\d*$'),  # Uppercase with numbers
            ],
            'priority': 1  # High priority for species-specific genes
        },
    }
    
    def __init__(self):
        self.gene_identifiers = []
        self.identifier_type = None
    
    def analyze_identifiers(self, identifiers: List[str]) -> Dict[str, float]:
        """
        Analyze a list of gene identifiers to infer species using progressive sampling.
        
        Args:
            identifiers: List of gene identifier strings
            
        Returns:
            Dictionary mapping species names to confidence scores
        """
        # For identifier type detection, use a sample to avoid bias
        if len(identifiers) > 100:
            sample_identifiers = random.sample(identifiers, 100)
        else:
            sample_identifiers = identifiers
        
        # Temporarily set sample for identifier type detection
        self.gene_identifiers = sample_identifiers
        self._determine_identifier_type()
        
        # Set all identifiers for progressive sampling
        self.gene_identifiers = identifiers
        
        # Progressive sampling for species inference
        return self._infer_species()
    
    def _determine_identifier_type(self):
        """Determine the type of gene identifiers (Ensembl, Entrez, or Symbol)."""
        if not self.gene_identifiers:
            raise SpeciesInferenceError("No gene identifiers provided")
        
        ensembl_count = 0
        entrez_count = 0
        symbol_count = 0
        
        for identifier in self.gene_identifiers:
            if self._is_ensembl_id(identifier):
                ensembl_count += 1
            elif self._is_entrez_id(identifier):
                entrez_count += 1
            else:
                symbol_count += 1
        
        # Determine the most common type
        if ensembl_count > entrez_count and ensembl_count > symbol_count:
            self.identifier_type = 'ensembl'
        elif entrez_count > symbol_count:
            self.identifier_type = 'entrez'
        else:
            self.identifier_type = 'symbol'
    
    def _is_ensembl_id(self, identifier: str) -> bool:
        """Check if an identifier matches any Ensembl pattern."""
        for pattern in self.ENSEMBL_PATTERNS.values():
            if pattern.match(identifier):
                return True
        return False
    
    def _is_entrez_id(self, identifier: str) -> bool:
        """Check if an identifier matches Entrez ID pattern."""
        return bool(re.match(r'^\d+$', identifier))
    
    def _infer_species(self) -> Dict[str, float]:
        """Infer species based on gene identifiers."""
        if self.identifier_type == 'ensembl':
            return self._infer_from_ensembl()
        elif self.identifier_type == 'entrez':
            return self._infer_from_entrez()
        else:
            return self._infer_from_symbols()
    
    def _infer_from_ensembl(self) -> Dict[str, float]:
        """Infer species from Ensembl IDs."""
        species_scores = {}
        
        for species, pattern in self.ENSEMBL_PATTERNS.items():
            matches = sum(1 for gene_id in self.gene_identifiers if pattern.match(gene_id))
            if matches > 0:
                species_scores[species] = matches / len(self.gene_identifiers)
        
        return species_scores
    
    def _infer_from_entrez(self) -> Dict[str, float]:
        """Infer species from Entrez IDs using range-based analysis."""
        species_scores = {}
        
        # Count numeric identifiers
        numeric_ids = [gene_id for gene_id in self.gene_identifiers if gene_id.isdigit()]
        
        if not numeric_ids:
            return species_scores
        
        # Convert to integers
        ids = [int(gene_id) for gene_id in numeric_ids]
        
        # Count how many IDs fall into each species range
        for species, (min_range, max_range) in self.ENTREZ_RANGES.items():
            matches = sum(1 for gene_id in ids if min_range <= gene_id <= max_range)
            if matches > 0:
                species_scores[species] = matches / len(numeric_ids)
        
        # If no IDs fall into expected ranges, use fallback logic
        if not species_scores:
            # Analyze the range of IDs for educated guesses
            min_id = min(ids)
            max_id = max(ids)
            
            if max_id < 100000:
                species_scores['human'] = 0.6
                species_scores['mouse'] = 0.3
                species_scores['rat'] = 0.1
            elif max_id < 200000:
                species_scores['mouse'] = 0.5
                species_scores['human'] = 0.3
                species_scores['rat'] = 0.2
            elif max_id < 300000:
                species_scores['rat'] = 0.5
                species_scores['mouse'] = 0.3
                species_scores['human'] = 0.2
            else:
                # Higher ranges might be other species
                species_scores['zebrafish'] = 0.4
                species_scores['fly'] = 0.3
                species_scores['worm'] = 0.3
        
        return species_scores
    
    def _infer_from_symbols(self) -> Dict[str, float]:
        """Infer species from gene symbols using optimized set operations."""
        species_scores = {}
        total_genes = len(self.gene_identifiers)
        
        # Convert gene identifiers to set for O(1) lookup
        gene_set = set(self.gene_identifiers)
        
        # Step 1: Check for exact matches using set intersection (much faster)
        exact_matches_found = False
        
        for species, data in self.SPECIES_GENES.items():
            symbols = data['symbols']
            # Use set intersection - O(min(n,m)) instead of O(n*m)
            exact_matches = len(gene_set & symbols)
            
            if exact_matches > 0:
                exact_matches_found = True
                # High confidence score for exact matches
                species_scores[species] = exact_matches * 10 / total_genes
        
        # Step 2: If no exact matches found, use pattern matching
        if not exact_matches_found:
            # Sort species by priority (most distinctive first)
            sorted_species = sorted(self.SPECIES_GENES.items(), key=lambda x: x[1]['priority'])
            
            for species, data in sorted_species:
                patterns = data['patterns']
                priority = data['priority']
                
                # Count pattern matches - iterate once through genes, check all patterns
                pattern_matches = 0
                for gene_id in self.gene_identifiers:
                    # Check if any pattern matches (short-circuit on first match)
                    if any(pattern.match(gene_id) for pattern in patterns):
                        pattern_matches += 1
                
                if pattern_matches > 0:
                    # Calculate score with priority weighting
                    priority_weight = (12 - priority) / 10
                    total_score = pattern_matches * priority_weight / total_genes
                    species_scores[species] = total_score
        
        # Step 3: If still no matches, provide fallback
        if not species_scores:
            # Default to human as most common species in research
            species_scores['human'] = 0.1  # Low confidence score
        
        return species_scores


class CountMatrixAnalyzer:
    """Analyzes count matrices to determine structure and extract gene identifiers."""
    
    def __init__(self, file_path: str):
        self.file_path = Path(file_path)
        self.df = None
        self.gene_orientation = None
        self.gene_identifiers = []
    
    def analyze(self) -> List[str]:
        """
        Analyze the count matrix and return gene identifiers.
        
        Returns:
            List of gene identifiers
        """
        self._load_data()
        self._determine_gene_orientation()
        self._extract_gene_identifiers()
        return self.gene_identifiers
    
    def _load_data(self):
        """Load the data file using polars."""
        if not self.file_path.exists():
            raise SpeciesInferenceError(f"File not found: {self.file_path}")
        
        # Determine file format
        if self.file_path.suffix.lower() == '.csv':
            self.df = pl.read_csv(self.file_path, has_header=True)
        elif self.file_path.suffix.lower() in ['.tsv', '.txt']:
            self.df = pl.read_csv(self.file_path, separator='\t', has_header=True)
        else:
            # Try CSV first, then TSV
            try:
                self.df = pl.read_csv(self.file_path, has_header=True)
            except Exception:
                try:
                    self.df = pl.read_csv(self.file_path, separator='\t', has_header=True)
                except Exception as e:
                    raise SpeciesInferenceError(f"Could not read file: {e}")
    
    def _determine_gene_orientation(self):
        """Determine if genes are in rows or columns."""
        if self.df is None:
            raise SpeciesInferenceError("Data not loaded")
        
        n_rows, n_cols = self.df.shape
        
        # Get first column name and first few values
        first_col = self.df.columns[0]
        first_col_values = self.df.select(pl.col(first_col)).head(10).to_series().to_list()
        
        # Check if first column contains gene identifiers vs cell barcodes
        # Gene identifiers: short names like "ACT1", "TP53", "YDL247W-A"
        # Cell barcodes: long strings like "AAACCTGTCGGAAACG-2"
        
        gene_like_count = 0
        cell_barcode_like_count = 0
        
        for val in first_col_values:
            if isinstance(val, str):
                # Check for cell barcode patterns (long strings with dashes and numbers)
                if len(val) > 15 and '-' in val and any(c.isdigit() for c in val):
                    cell_barcode_like_count += 1
                # Check for generic cell identifiers (cell1, cell2, etc.)
                elif val.startswith('cell') and val[4:].isdigit():
                    cell_barcode_like_count += 1
                # Check for gene identifier patterns (shorter, more gene-like)
                elif len(val) <= 15 and not val.replace('.', '').replace('-', '').isdigit():
                    gene_like_count += 1
        
        # If we have more cell barcode-like patterns, genes are likely in columns
        if cell_barcode_like_count > gene_like_count:
            self.gene_orientation = 'columns'
        else:
            # Check if first row contains gene identifiers
            first_row = self.df.head(1).to_pandas().iloc[0].tolist()
            gene_like_count_row = sum(1 for val in first_row[1:]  # Skip first element
                                    if isinstance(val, str) and len(val) <= 15 and not val.replace('.', '').replace('-', '').isdigit())
            
            if gene_like_count_row >= len(first_row) * 0.7:
                self.gene_orientation = 'columns'
            else:
                # Default to rows if uncertain
                self.gene_orientation = 'rows'
    
    def _extract_gene_identifiers(self):
        """Extract gene identifiers based on orientation."""
        if self.gene_orientation == 'rows':
            # Genes are in rows, first column contains identifiers
            self.gene_identifiers = self.df.select(pl.col(self.df.columns[0])).to_series().to_list()
        else:
            # Genes are in columns, column names contain identifiers
            # Skip the first column (usually barcode or index)
            self.gene_identifiers = self.df.columns[1:]
        
        # Filter out None values and convert to strings
        self.gene_identifiers = [str(gene_id) for gene_id in self.gene_identifiers if gene_id is not None]


def infer_species_from_mtx(features_file: str) -> Tuple[str, int]:
    """
    Infer species from a 10X Genomics features file (TSV format).
    
    Args:
        features_file: Path to the features.tsv.gz file
        
    Returns:
        Tuple of (inferred species name, number of genes)
    """
    try:
        # Load the features file
        features_path = Path(features_file)
        
        # Determine file format and read accordingly
        if features_path.suffix.lower() == '.gz':
            # Compressed file
            if features_path.stem.endswith('.tsv'):
                df = pl.read_csv(features_file, separator='\t', has_header=False)
            else:
                df = pl.read_csv(features_file, has_header=False)
        else:
            # Uncompressed file
            if features_path.suffix.lower() in ['.tsv', '.txt']:
                df = pl.read_csv(features_file, separator='\t', has_header=False)
            else:
                df = pl.read_csv(features_file, has_header=False)
        
        # Get number of genes (rows) from features file
        num_genes = df.height
        
        # Extract gene identifiers from first column (Ensembl IDs)
        gene_identifiers = df.select(pl.col(df.columns[0])).to_series().to_list()
        
        if not gene_identifiers:
            raise SpeciesInferenceError("No gene identifiers found in the features file")
        
        # Analyze gene identifiers to infer species
        gene_analyzer = GeneIdentifierAnalyzer()
        species_scores = gene_analyzer.analyze_identifiers(gene_identifiers)
        
        if not species_scores:
            raise SpeciesInferenceError("Could not infer species from gene identifiers")
        
        # Get the best species using tie-breaking logic
        best_species = _select_best_species(species_scores)
        
        # Convert to standard format (e.g., "homo-sapiens")
        species = SPECIES_NAME_MAPPING.get(best_species, best_species)
        return (species, num_genes)
        
    except Exception as e:
        raise SpeciesInferenceError(f"Error inferring species from MTX features file: {e}")


def _infer_species_from_anndata(adata) -> str:
    """
    Common function to infer species from an AnnData object.
    
    Args:
        adata: AnnData object
        file_type: String describing the file type (for error messages)
        
    Returns:
        Inferred species name
    """
    try:
        # Extract gene identifiers from var_names (genes)
        gene_identifiers = list(adata.var_names)
        
        if not gene_identifiers:
            raise SpeciesInferenceError(f"No gene identifiers found")
        
        # Analyze gene identifiers to infer species
        gene_analyzer = GeneIdentifierAnalyzer()
        species_scores = gene_analyzer.analyze_identifiers(gene_identifiers)
        
        if not species_scores:
            raise SpeciesInferenceError("Could not infer species from gene identifiers")
        
        # Get the best species using tie-breaking logic
        best_species = _select_best_species(species_scores)
        
        # Convert to standard format
        return SPECIES_NAME_MAPPING.get(best_species, best_species)
        
    except Exception as e:
        raise SpeciesInferenceError(f"Error inferring species: {e}")

def infer_species(input_file: str) -> str:
    """
    Main function to infer species from a count matrix file.
    
    Args:
        input_file: Path to the input CSV/TSV file
        
    Returns:
        Inferred species name
    """
    try:
        # Analyze the count matrix
        matrix_analyzer = CountMatrixAnalyzer(input_file)
        gene_identifiers = matrix_analyzer.analyze()
        
        if not gene_identifiers:
            raise SpeciesInferenceError("No gene identifiers found in the file")
        
        # Analyze gene identifiers to infer species
        gene_analyzer = GeneIdentifierAnalyzer()
        species_scores = gene_analyzer.analyze_identifiers(gene_identifiers)
        
        if not species_scores:
            raise SpeciesInferenceError("Could not infer species from gene identifiers")
        
        # Get the best species using tie-breaking logic
        best_species = _select_best_species(species_scores)
        
        # Convert to standard format (e.g., "homo-sapiens")
        return SPECIES_NAME_MAPPING.get(best_species, best_species)
        
    except Exception as e:
        raise SpeciesInferenceError(f"Error inferring species: {e}")


def main():
    """Main function to handle command line arguments and execute the analysis."""
    parser = argparse.ArgumentParser(
        description="Infer species from single-cell RNA-seq count matrix",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Examples:
    python infer_species.py data.csv --format csv
    python infer_species.py data.tsv --format csv
    python infer_species.py features.tsv.gz --format mtx
    python infer_species.py data.h5ad --format h5ad
    python infer_species.py data.h5ad --format h5ad-multi-sample
    python infer_species.py data.h5 --format h5
        """
    )
    
    parser.add_argument(
        'input_file',
        help='Path to the input file (CSV/TSV for csv format, TSV for mtx format, h5ad for h5ad format, h5 for h5 format)'
    )
    
    parser.add_argument(
        '--format', '-f',
        choices=['csv', 'mtx', 'h5ad', 'h5ad-multi-sample', 'h5'],
        default='csv',
        help='Input format: csv for count matrix files, mtx for 10X Genomics features file, h5ad for AnnData files, h5ad-multi-sample for multi-sample AnnData files, h5 for 10X Genomics Cell Ranger HDF5 files (default: csv)'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='species.txt',
        help='Output file name (default: species.txt)'
    )
    
    parser.add_argument(
        '--barcodes',
        help='Path to barcodes.tsv file (required for mtx format to calculate cell count)'
    )
    
    args = parser.parse_args()
    
    try:
        # Infer species based on format
        if args.format == 'csv':
            species = infer_species(args.input_file)
            
            # Detect and write gene format
            analyzer = CountMatrixAnalyzer(args.input_file)
            gene_identifiers = analyzer.analyze()  # This initializes the dataframe
            gene_analyzer = GeneIdentifierAnalyzer()
            gene_analyzer.analyze_identifiers(gene_identifiers)
            
            # Map identifier type to format name
            gene_format = GENE_FORMAT_MAPPING.get(gene_analyzer.identifier_type, 'unknown')
            
            # Calculate data size: rows = genes, columns = cells
            n_rows, n_cols = analyzer.df.shape
            if analyzer.gene_orientation == 'rows':
                # Genes in rows, cells in columns (excluding first column which is gene names)
                num_genes = n_rows
                num_cells = n_cols - 1  # Subtract the gene identifier column
            else:
                # Genes in columns, cells in rows (excluding first column which is cell names)
                num_genes = n_cols - 1  # Subtract the cell identifier column
                num_cells = n_rows
            
        elif args.format == 'mtx':
            species, num_genes = infer_species_from_mtx(args.input_file)
            
            # For MTX format, we know it's Ensembl IDs from the features file
            gene_format = 'Ensembl Id'
            
            # Number of cells (columns) from barcodes file
            if args.barcodes:
                barcodes_path = Path(args.barcodes)
                if barcodes_path.suffix.lower() == '.gz':
                    if barcodes_path.stem.endswith('.tsv'):
                        barcodes_df = pl.read_csv(args.barcodes, separator='\t', has_header=False)
                    else:
                        barcodes_df = pl.read_csv(args.barcodes, has_header=False)
                else:
                    if barcodes_path.suffix.lower() in ['.tsv', '.txt']:
                        barcodes_df = pl.read_csv(args.barcodes, separator='\t', has_header=False)
                    else:
                        barcodes_df = pl.read_csv(args.barcodes, has_header=False)
                num_cells = barcodes_df.height
            else:
                # Number of cells (columns) is unknown from features file alone
                num_cells = 0
            
        elif args.format in ['h5ad', 'h5ad-multi-sample', 'h5']:
            if args.format == 'h5':
                adata = sc.read_10x_h5(args.input_file)
            else:
                adata = sc.read_h5ad(args.input_file)
            
            # Infer species using the already-loaded adata
            species = _infer_species_from_anndata(adata)
            
            # Detect gene format from the same adata object
            gene_identifiers = list(adata.var_names)
            
            gene_analyzer = GeneIdentifierAnalyzer()
            gene_analyzer.analyze_identifiers(gene_identifiers)
            
            gene_format = GENE_FORMAT_MAPPING.get(gene_analyzer.identifier_type, 'unknown')
            
            # Calculate data size: adata.shape = (n_obs, n_vars) = (cells, genes)
            # So rows = genes, columns = cells
            num_cells, num_genes = adata.shape
        
        # Write species output
        output_path = Path(args.output)
        with open(output_path, 'w') as f:
            f.write(species)
        
        # Write gene format output
        format_output_path = Path('gene_format.txt')
        with open(format_output_path, 'w') as f:
            f.write(gene_format)
        
        # Write data size output (rows columns)
        data_size_output_path = Path('data_size.txt')
        with open(data_size_output_path, 'w') as f:
            f.write(f"{num_genes} {num_cells}")
        
        print(f"Species inferred: {species}")
        print(f"Gene format detected: {gene_format}")
        print(f"Data size: {num_genes} genes (rows) × {num_cells} cells (columns)")
        print(f"Species result written to: {output_path}")
        print(f"Gene format result written to: {format_output_path}")
        print(f"Data size result written to: {data_size_output_path}")
        
    except SpeciesInferenceError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
