#!/usr/bin/python
# -*- coding: UTF-8 -*-
# search for GAACCUGAACCU; AACCUGAACCUG; ACCUGAACCUGA; CCUGAACCUGAA; CUGAACCUGAAC; UGAACCUGAACC. 
import sys
def fun(argv):
    input_file = argv[0]
    tagged_file = argv[1]
    untagged_file = argv[2]
    print 'input file:   %s' % (input_file)
    print 'search range: fist 50nt'
    print 'subsequence:  GAACCUGAACCU(12nt)' 
    print 'tagged file:  %s' % (tagged_file)
    print 'untagged file:%s' % (untagged_file)
    records = []
    with open(input_file, 'r') as f:
        record = []
        for line in f.readlines():
            record.append(line.strip())
            if len(record) == 4:
                records.append(record)
                record = []
    tagged_count = 0
    untagged_count = 0
    with open(tagged_file, 'w') as f, open (untagged_file, 'w') as g:
        for record in records:
            sequence = record[1][:50]
            if sequence.find('GAACCUGAACCU') < 0 and sequence.find('AACCUGAACCUG') < 0 and sequence.find('ACCUGAACCUGA') < 0 and sequence.find('CCUGAACCUGAA') < 0 and sequence.find('CUGAACCUGAAC') < 0 and sequence.find('UGAACCUGAACC') < 0:
            	untagged_count = untagged_count + 1
		# print untagged sequence
	        for r in record:
		    if r == '':
		        continue
		    g.write(r + '\n')
		continue
            # print tagged sequence
            tagged_count = tagged_count + 1
            for r in record:
                if r == '':
                    continue
                f.write(r + '\n')
    print 'total [%d] tagged [%d] untagged [%d]' % (len(records), tagged_count, untagged_count)
if __name__ == '__main__':
    run_format  = '  python main.py input.file tagged.file untagged.file'
    run_example = '  python main.py input.fastq tagged.file untagged.file'
    print len(sys.argv)
    if len(sys.argv) != 4:
        print 'you need follow this run format'
        print run_format
        print run_example
        exit(1)
    sys.argv
    fun(sys.argv[1:])
