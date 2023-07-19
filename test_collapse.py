def collapse_reads(infile, outfile):
        counts = []
        reads = []
        
        with open(stdoutfile, 'r') as f, open(output, 'w') as g:
            for line in f:
               if not line.startswith('>') and not line in reads:
                   reads.append(line)
                   counts.append(1)
               elif not line.startswith('>') and line in reads:
                   counts(reads.index(line)) += 1

            d = OrderedDict(zip(reads, counts))
            reads_tmp = list(OrderedDict(sorted(d.items(), key=lambda x: x[1])))
            reads_final = reads_tmp[::-1]
            counts.sort(reverse = True)
            
            for idx, line in enumerate(reads_final):
                g.write('>Read-' + str(counts[idx]))
                
                if not umi:
                    g.write(line)
                elif umi > 0:
                    g.write(line[umi-1:])
                elif umi < 0:
                    g.write(line[0:umi])
        f.close()
        g.close()
