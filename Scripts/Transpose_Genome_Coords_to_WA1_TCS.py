import argparse
import re
cigar_regex = re.compile(r"[^\W\d_]+|\d+")
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("Changes", help= "Input BED File")
parser.add_argument("InBED", help= "Input BED File")
parser.add_argument("OutBED", help= "Root of output files")
args = parser.parse_args()
## Handle options

Dict = {}
with open(str(args.Changes),'r') as In:
    data = In.readline().rstrip()
    while data:
        data = cigar_regex.findall(data)
 #       print(data)
        if len(data[0]) > 1:
            #Deletion
            Coord = int(data[1])
            DelSize = len(data[0]) - 1
            Dict[Coord] = DelSize
        if len(data[2]) > 1:
            #Insertion
            Coord = int(data[1])
            InsSize = len(data[2]) - 1
            Dict[Coord] = -InsSize
        data = In.readline().rstrip()
        
#print(Dict)

Out = open(str(args.OutBED),'w')
with open(str(args.InBED),'r') as In:
    Header = In.readline()
    Out.write(Header)
    data = In.readline().rstrip().split()
    while data:
            StartCoord = int(data[1])
            StopCoord = int(data[2])
            NewStardCoord = StartCoord
            NewStopCoord = StopCoord
            for i in Dict:
                if StartCoord > i:
                    NewStardCoord += Dict[i]
                if StopCoord > i:
                    NewStopCoord += Dict[i]
  #          print(StartCoord, NewStardCoord, StopCoord, NewStopCoord)
            NewLine = [data[0]] + [str(NewStardCoord)] + [str(NewStopCoord)] + data[3:]
            Out.write('\t'.join(NewLine) + '\n')
            data = In.readline().rstrip().split()
Out.close()

            
                
            
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                

        
    
    