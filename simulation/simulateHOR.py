import random
import os
def buildHOR_monomer_seqs(hor_monomer_len,hor_number,type):
    hors_monomer_seq = []
    if type[0] == 'canonical':
        for i in range(hor_number):
            hor_monomer_seq = []
            for j in range(hor_monomer_len):
                hor_monomer_seq.append(j + 1)
            hors_monomer_seq.append(hor_monomer_seq)
    elif type[0] == 'LN':
        LN_number = type[1]
        random_LN_hors_monomer_seq = []
        for i in range(LN_number):
            hor_monomer_seq = []
            for j in range(hor_monomer_len):
                hor_monomer_seq.append(j + 1)
            LN_unit_len = random.randint(1,4)
            LN_unit_copynumber = random.randint(2,5)
            LN_location = random.randint(0,hor_monomer_len - LN_unit_len)
            prefix = hor_monomer_seq[0:LN_location]
            surfix = hor_monomer_seq[LN_location + LN_unit_len:]
            LN_mono_seq = []
            LN_mono_seq_unit = hor_monomer_seq[LN_location:LN_location+LN_unit_len]
            for j in range(LN_unit_copynumber):
                LN_mono_seq += LN_mono_seq_unit
            random_LN_hors_monomer_seq.append(prefix + LN_mono_seq + surfix)
        for i in range(hor_number - LN_number):
            hor_monomer_seq = []
            for j in range(hor_monomer_len):
                hor_monomer_seq.append(j + 1)
            random_LN_hors_monomer_seq.append(hor_monomer_seq)
        random.shuffle(random_LN_hors_monomer_seq)
        hors_monomer_seq = random_LN_hors_monomer_seq
    return hors_monomer_seq

based = ['A','T','C','G']
template_seq_len_list = [100,200,300,400]
monomer_mutation_rate_list = [0.3,0.2,0.1]
Hor_mutation_rate_list = [0.025,0.015,0.005]

hor_unit_len = 5
hor_number = 40
canonical_type = ('canonical',0)
LN_type = ('LN',20)

workdir = './canonical'
# random monomer sequence
hors_monomer_seq = buildHOR_monomer_seqs(hor_unit_len,hor_number,canonical_type)
outmonomer_seq_file = workdir + '/monomer_sequence.txt'
outmonomer_seq_file = open(outmonomer_seq_file,'w')
print(hors_monomer_seq)
for hor in hors_monomer_seq:
    hor_DNA_seq = []
    for mo in hor:
        outmonomer_seq_file.write(str(mo)+' ')
outmonomer_seq_file.close()

# random DNA sequence
for i in template_seq_len_list:
    for j in monomer_mutation_rate_list:
        for k in Hor_mutation_rate_list:
            template_seq_len = i
            monomer_mutation_rate = j
            Hor_mutation_rate = k
            print(str(template_seq_len)+'\t'+str(monomer_mutation_rate)+'\t'+str(Hor_mutation_rate))
            outdir = workdir + '/' + 'test.' + str(template_seq_len)+'.'+\
                                   str(monomer_mutation_rate)+'.'+\
                                   str(Hor_mutation_rate)
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            template_seq = []
            for b in range(i):
                random.shuffle(based)
                template_seq.append(based[0])
            outtemplate_seq_file = outdir + '/template_seq.fa'
            outtemplate_seq_file = open(outtemplate_seq_file, 'w')
            outtemplate_seq_file.write('>template'+str(template_seq_len)+'.'+\
                                   str(monomer_mutation_rate)+'.'+\
                                   str(Hor_mutation_rate) +'\n')
            for b in template_seq:
                outtemplate_seq_file.write(str(b))
            outtemplate_seq_file.write('\n')
            outtemplate_seq_file.close()
            hor_unit_seqs = {}
            for mo in range(hor_unit_len):
                mutation_index_list = list(range(template_seq_len))
                random.shuffle(mutation_index_list)
                mutation_index = mutation_index_list[:int(template_seq_len*monomer_mutation_rate)]
                monomer_seq = template_seq.copy()
                for m in mutation_index:
                    new_based = []
                    for b in based:
                        if b != monomer_seq[m]:
                            new_based.append(b)
                    random.shuffle(new_based)
                    monomer_seq[m] = new_based[0]
                hor_unit_seqs[mo+1] = monomer_seq

            for mo in hor_unit_seqs.keys():
                outmonomer_file = outdir + '/monomer_seq'+str(mo)+ '.fa'
                outmonomer_file = open(outmonomer_file,'w')
                outmonomer_file.write('>'+str(mo)+'\n')

                for b in hor_unit_seqs[mo]:
                    outmonomer_file.write(str(b))
                outmonomer_file.write('\n')
                outmonomer_file.close()
            HORs_DNA_seq = []
            for hor in hors_monomer_seq:
                hor_DNA_seq = []
                for mo in hor:
                    mutation_index_list = list(range(template_seq_len))
                    random.shuffle(mutation_index_list)
                    mutation_index = mutation_index_list[:int(template_seq_len * Hor_mutation_rate)]
                    monomer_seq = hor_unit_seqs[mo].copy()
                    for m in mutation_index:
                        new_based = []
                        for b in based:
                            if b != monomer_seq[m]:
                                new_based.append(b)
                        random.shuffle(new_based)
                        monomer_seq[m] = new_based[0]
                    hor_DNA_seq.append(monomer_seq)
                HORs_DNA_seq.append(hor_DNA_seq)
            outHOR_fa_file = outdir + '/HOR.fa'
            outHOR_fa_file = open(outHOR_fa_file,'w')
            outHOR_fa_file.write('>HOR\n')
            for hor in HORs_DNA_seq:
                for mo in hor:
                    for b in mo:
                        outHOR_fa_file.write(b)
            outHOR_fa_file.write('\n')
            outHOR_fa_file.close()

