# Lisa Boatner
# Date Created: 11/19/20
# Date Modified: 230314
# Notes: Universal AA automated "all" parameter

import os, sys, csv, re, argparse, re, collections, glob, statistics, numbers, string, time
from os import listdir, path
from os.path import isfile, isdir, join
from Bio import SeqIO
import pandas as pd

# check a file exists in current directory
def check_file_exists(file):
  if (path.exists(file) == False):
    print(file + " does not exist in " + str(os.getcwd()))
    sys.exit()

def read_fasta(file, anon_dict):
  check_file_exists(file)
  for record in SeqIO.parse(file, "fasta"):
    anon_dict[record.id] = record.seq
  return anon_dict

def get_aa_identifiers(sequence, protein, aa):
  temp_lst = []
  pro_lst = []
  for i in range(len(sequence)):
    if sequence[i] == aa:
      temp_lst.append(protein + '_' + aa + str(i + 1))
      pro_lst.append(protein)

  return [temp_lst, pro_lst]

def list_to_string(lst):
  return (';'.join([str(elem) for elem in lst]))

def write_file(file, header, misc, special, args):
  out_file = open(file, 'w')

  header = 'protein,' + args.aa.lower() + '_ids'
  out_file.write(header + '\n')
  
  if isinstance(misc, dict):
    write_dict_file(out_file, misc, special)

  elif isinstance(misc, list):
    write_list_file(out_file, misc, special)

  else:
    print("Cannot write file, unrecognized format. " + str(misc))

  out_file.close()

def write_dict_file(out_file, misc, special):
  for k in misc:
    peptides = misc[k]
    for i in range(len(peptides)):
      out_file.write(str(k) + ',' + str(peptides[i]) + '\n')

def write_list_file(out_file, misc, special):
  for i in range(len(misc)):
    current = misc[i]
    if special == True:
      st = current.output()
      out_file.write(st + '\n')
    else:
      st = ''
      for j in range(len(current)):
        st +=  str(current[j]) + ','
      out_file.write(st[:-1] + '\n')

def start_game(args, cd, aad):
  uniprot_dict = {}
  uniprot_dict = read_fasta(args.i, uniprot_dict)

  count = 0
  identifiers_dict = {}

  proteins = []
  aas = []

  for k in uniprot_dict:
    protein_seq = uniprot_dict[k]
    protein = k.split('|')[1]
    results = get_aa_identifiers(protein_seq, protein, args.aa)
    proteins += results[1]
    aas += results[0]

  if args.wo == "True":
    # proteinid, cysteineid
    new_df = pd.DataFrame()
    new_df['proteinid'] = proteins
    new_df['residueid'] = aas
    os.chdir(aad)
    new_df.to_csv(args.o, index = False)
    os.chdir(cd)

def auto_start(args, aa_dict):

  cd = os.getcwd()

  path_data = os.path.join(os.getcwd(), 'aa_ids')
  if not os.path.exists(path_data):
    os.makedirs(path_data)

  if args.all == 'True':
    print("running all amino acids")
    aas = list(aa_dict.keys())
    resnames = list(aa_dict.values())

    for i in range(len(aas)):
      print("running " + aas[i])
      args.aa = aas[i]
      args.rid = resnames[aas.index(aas[i])]
      args.o = '2401_uniprot_' + args.rid + 'ids.csv'
      start_game(args, cd, path_data)
  else:
    print("running " + args.rid)
    args.o = '2401_uniprot_' + args.rid + 'ids.csv'
    start_game(args, cd, path_data)

def main():

  aa_dict = {'A': 'alanine', 'G': 'glycine', 'I': 'isoleucine', 'L': 'leucine', 'P': 'proline', 'V': 'valine',
  'F': 'phenylalanine', 'W': 'tryptophan', 'Y': 'tyrosine', 'D': 'aspartate', 'E': 'glutamate', 'R': 'arginine',
  'H': 'histidine', 'K': 'lysine', 'S': 'serine', 'T': 'threonine', 'C': 'cysteine', 'M': 'methionine', 'N': 'asparagine',
  'Q': 'glutamine'}

  # https://www.dropbox.com/s/r8miefs6m1sb2it/2301_uniprot.fasta?dl=1

  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input', dest='i', nargs='?', default="2401_uniprot.fasta", type=str, help='default 2401_uniprot.fasta')
  parser.add_argument('-all', '--all', dest='all', nargs='?', default="False", type=str, help='default False, options: True')
  parser.add_argument('-aa', '--amino_acid', dest='aa', nargs='?', default="C", type=str, help='default C')
  parser.add_argument('-rid', '--resname', dest='rid', nargs='?', default="cysteine", type=str, help='default cysteine')
  parser.add_argument('-o', '--output', dest='o', nargs='?', default="2401_uniprot_ids.csv", type=str, help='default 2401_uniprot_ids.csv')
  parser.add_argument('-wo', '--write_output', dest='wo', nargs='?', default="True", type=str, help='default True')
    
  args = parser.parse_args()
  auto_start(args, aa_dict)

if __name__ == "__main__":
  main()
