from Bio import SeqIO

a2m_file = str(snakemake.input)
out_file = str(snakemake.output)

with open(a2m_file) as fh:
    a2m = SeqIO.parse(fh, 'fasta')
    ss_pred = next(a2m)
    ss_conf = next(a2m)
    cons = next(a2m)

    num_pos = len(cons)
    gap_nums = [0] * num_pos

    num_records = 0
    for record in a2m:
        for i in range(num_pos):
            gap_nums[i] += record.seq[i] == '-'
        num_records += 1

with open(out_file, 'w') as out:
    for i in range(num_pos):
        gap_score = gap_nums[i] / num_records
        sec_struct = ss_pred.seq[i]
        if sec_struct == 'H' and gap_score <= 0.1:
            out.write(str(i+1) + "\n")
