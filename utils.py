import requests
import pandas as pd
from Bio.Blast import NCBIWWW
import re


def find_protein_seq(essential_genes):
    """Takes list of  genes and returns their protein sequence, using BIGG database"""
    protein_seqs = {}
    for gene in essential_genes:
        #print(gene)
        try:
            text = requests.get(f'http://bigg.ucsd.edu/models/iIT341/genes/{gene}').text
            ls = text.split('\n')
            string = ''
            for i in range(len(ls)):
                if '<h4>Protein Sequence</h4>' in ls[i]:
                    string = ls[i + 1]
                    break
            protein_seqs[gene] = string.replace('          <p class="sequence">', '').replace('</p>', '')
        except:
            protein_seqs[gene] = None
    return protein_seqs


def write_fasta(protein_seqs, filename):
    """Writes fasta file based on protein sequences."""
    string = ''
    for gene in protein_seqs:
        string += f'>{gene}\n{protein_seqs[gene]}\n'
    file = open(filename, 'w')
    file.write(string)
    file.close()


def read_table_genes_to_validate(filename):
    """Reads csv file and returns all experimentally validated essential genes. """
    df = pd.read_csv(filename)
    ls = df.query("essentiality=='E'")["locus"]
    
    return ls

def build_dictionary(filename):
    """Builds a dictionary with the IDs and their e-values"""
    file = open(filename, 'r')
    lines = file.read().split('\n')
    file.close()
    dictionary = {}
    key = ''
    for line in lines:
        if '  <Iteration_query-def>' in line:
            key = line.replace('  <Iteration_query-def>', '').replace('</Iteration_query-def>', '')
            dictionary[key] = []
        elif '      <Hsp_evalue>' in line:
            dictionary[key].append(float(line.replace('      <Hsp_evalue>', '').replace('</Hsp_evalue>', '')))
    return dictionary


def find_non_homologous(filename, evalue):
    """Analyses a dictionary and return a list of IDs with higher or equal of the e-value or with no hits"""
    dictionary = build_dictionary(filename)
    ls = []
    for key in dictionary:
        if len(dictionary[key]) == 0:
            ls.append(key)
        elif min(dictionary[key]) >= evalue:
            ls.append(key)
    return ls

#write_fasta(find_protein_seq(['HP0370', 'HP0950']))


def blast(file_name):
    # Perform BLASTp search for the sequence
    with open(file_name,'r') as f:
        for line in f:
            if line.startswith('>'):
                query_id = line.strip()[1:]
                query_seq = next(f).strip()
                results = NCBIWWW.qblast('blastp', 'nr', query_seq, entrez_query='Homo sapiens [organism]', expect= 0.001, hitlist_size=3)
                result_file = 'result_file.xml'
                with open(result_file, "w") as result_f:
                    result_f.write(results.read())
                print(f"BLASTP results for query {query_id} are saved in {result_file}")
                result_f.close()
    f.close()

#blast("FastaEssential.fasta")


def retrieve_id(file_name):
    hits = 0
    id = []
    with open(file_name, 'r') as f:
        for line in f:
            if line.startswith('>'):
                hits +=1
                query_id = line.strip()[1:]
                query_seq = next(f).strip()
                id.append(query_id)
    return id

def hits_drug_bank(file_name):
    #para ver os resultados do drug bank
    pattern_results = '\([1-9][0-9]*\s+result(s)?\)|\([1-9][0-9]*\s+results\)'
    prot_id = []
    with open(file_name, 'r') as f:
        for line in f:
            match = re.search(pattern_results, line)
            if match:
                prot_id.append(line[:6])
    return prot_id

