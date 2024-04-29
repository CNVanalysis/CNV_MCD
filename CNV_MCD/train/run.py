from train.main import detect_cnv
if __name__ == '__main__':
    chr_num = 21;
    for i in range(1, 51):

        num = str(i)
        train_bam_path = '/.. /sim' + num + '_.sort.bam'
        output_path = '/.. /sim' + num + "_result.txt"
        detect_cnv(train_bam_path, output_path, chr_num)
