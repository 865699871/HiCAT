import argparse
import os
def main():
    parser = argparse.ArgumentParser(description="HiCAT: automated annotation centromere")
    parser.add_argument("-i", "--input_fasta",help="centromere DNA sequence in fasta format, required",required=True)
    parser.add_argument("-t", "--monomer_template", help="monomer template DNA sequence in fasta format for stringdecomposer to build block, required",required=True)
    parser.add_argument("-o", "--output_dir", help="HiCAT output path default is ./HiCAT_out",default='./HiCAT_out',required=False)

    parser.add_argument("-ms", "--min_similarity", help="The lower bound for similarity threshold which used to remove edges in block graph, default is 0.94", type=float, default=0.94,required=False)
    parser.add_argument("-st", "--step", help="The similarity threshold iteratively increases from min_similarity to nearly 1 with a specific step, default is 0.005", type=float, default=0.005,required=False)
    parser.add_argument("-mh", "--max_hor_len", help="An upper bound for the length of the tandem repeat unit by default 40 monomers for improving efficiency", type=int, default=40,required=False)

    parser.add_argument("-sp", "--show_hor_number",help="Default visualized the top five HORs", type=int, default=5,required=False)
    parser.add_argument("-sn", "--show_hor_min_repeat_number",help="Default visualized the HORs with repeat numbers greater than 10", type=int, default=10,required=False)

    parser.add_argument("-th", "--thread",help="The number of threads, default is 1", type=int, default=1,required=False)

    args = parser.parse_args()

    input_fasta = args.input_fasta
    output_dir = args.output_dir
    monomer_template = args.monomer_template


    min_similarity = args.min_similarity
    step = args.step
    max_hor_len = args.max_hor_len

    show_hor_number = args.show_hor_number
    show_hor_min_repeat_number = args.show_hor_min_repeat_number
    thread = args.thread

    script_path = os.path.split(os.path.realpath(__file__))[0]

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    monomer_template_base_sequence = ''
    with open(monomer_template, 'r') as f:
        f.readline()
        monomer_template_base_sequence = f.readline()[:-1]
    overlap_threshold = 2 * len(monomer_template_base_sequence)
    if overlap_threshold > 500:
        # run SD
        cmd = 'python ' + \
              script_path+ '/stringdecomposer/bin/stringdecomposer' + ' ' + \
              input_fasta + ' ' + monomer_template + \
              ' -v ' + str(overlap_threshold) + ' -o ' + output_dir
    else:
        cmd = 'python ' + \
              script_path + '/stringdecomposer/bin/stringdecomposer' + ' ' + \
              input_fasta + ' ' + monomer_template + \
              ' -o ' + output_dir
    print('Run stringdecomposer\n')
    print(cmd)
    os.system(cmd)

    input_seq = ''
    header = ''
    with open(input_fasta,'r') as f:
        while True:
            line = f.readline()[:-1]
            if not line:
                break
            if line.startswith('>'):
                header = line
            else:
                input_seq += line
    out_input_file = output_dir + '/input_fasta.1.fa'
    out_input_file = open(out_input_file,'w')
    out_input_file.write(header+'\n')
    out_input_file.write(input_seq+'\n')
    out_input_file.close()

    # Run HiCAT HOR
    cmd = 'python ' + \
          script_path+'/HiCAT_HOR.py' + ' ' + \
          '-d' + ' ' + output_dir + '/final_decomposition.tsv' + ' ' + \
          '-b' + ' ' + output_dir + '/input_fasta.1.fa' + ' ' + \
          '-o' + ' ' + output_dir + ' ' + \
          '-s' + ' ' + str(min_similarity) + ' ' + \
          '-st' + ' ' + str(step) + ' ' + \
          '-m' + ' ' + str(max_hor_len) + ' ' + \
          '-sp' + ' ' + str(show_hor_number) + ' ' + \
          '-sn' + ' ' + str(show_hor_min_repeat_number) + ' ' + \
          '-t' + ' ' + str(thread)
    print('Run HiCAT HOR\n')
    print(cmd)
    os.system(cmd)










if __name__ == '__main__':
    main()