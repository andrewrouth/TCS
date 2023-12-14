#!/bin/python3
import numpy as np
import argparse
from textwrap import wrap
parser = argparse.ArgumentParser()
parser.add_argument("InputGenome", help="Input Ref")
parser.add_argument("Coverage", help="Coverage bedgraph")
parser.add_argument("OutputGenome", help="Desired name of padded Ref")
parser.add_argument("--MetaDataFile", help="Report Name")
#parser.add_argument("--Fraction", help="Report Name")
parser.add_argument("--MinCount", help="Default is 10")
args = parser.parse_args()
InputGenome = str(args.InputGenome)
Coverage = str(args.Coverage)
OutputGenome = str(args.OutputGenome)
if args.MinCount:
    MinCount = str(args.MinCount)
else:
    MinCount = 10
    
AnnoDict = {}  ##Keys = Root  ##Values = Sample
if args.MetaDataFile:
     with open(str(args.MetaDataFile), 'r') as In:
         Data = In.readline().rstrip()
         while Data:
             Root, Name = Data.split()
             AnnoDict[Root] = Name
             Data = In.readline().rstrip()
else:
    pass

if args.MinCount:
    MinCount = str(args.MinCount)
else:
    MinCount = 10

##Read in REF
OriginalReference = ''
with open(InputGenome, 'r') as FASTAIN:
    Header = FASTAIN.readline().rstrip()
#    print("Original Reference = ", Header)
    data = FASTAIN.readline().rstrip()
    while data:
        OriginalReference += data
        data = FASTAIN.readline().rstrip()
        
Output = np.array([i for i in OriginalReference])       
with open(Coverage, 'r') as In:
    data = In.readline()
    while data:
        data = data.split()
        Start = int(data[1])
        Stop = int(data[2])
        Cov = float(data[3])
        if Cov < MinCount:
            Output[Start:Stop] = 'N'
        else:
            pass
        data = In.readline()
        
Output = ''.join(Output)
Ns = Output.count('N')
Completeness = 1 - (Ns/len(Output))
Completeness = round((Completeness*100),1)
Output = wrap(Output,60)
Output = '\n'.join(Output)

if args.MetaDataFile:
    Root = InputGenome.split('.fasta')[0]
    Dir, Root = Root.split("/")
    OutputGenome = AnnoDict[Root] + '.fasta'
    Header = '>SARS-CoV-2/human/USA/' + AnnoDict[Root] + '/2021' + '\tCompleteness: ' + str(Completeness) + '%'
else:
    pass

Out = open(OutputGenome, 'w')
Out.write(Header + '\n' + Output + '\n')