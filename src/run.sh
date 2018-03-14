#!/bin/bash
for number in {0..0}
do
for sample in {0..0}
do
echo $number
#python util.py -d 0222 -m random -n 20000 -t $number -v -0.05
#python util.py -d 0222 -m random -n 20000 -t $number -v 0
#python util.py -d 0222 -m random -n 20000 -t $number -v 0.05
#python util.py -d 0222 -m random -n 20000 -t $number -v 0.10
#python util.py -d 0222 -m random -n 20000 -t $number -v 0.15
#python util.py -d 0222 -m random -n 20000 -t $number -v 0.20
#python util.py -d 0222 -m random -n 20000 -t $number -v 0.25
#python util.py -d 0222 -m random -n 20000 -t $number -v 0.30
python util.py -d 0222 -m random -n 20000 -t $number -v 0.35
python util.py -d 0222 -m random -n 20000 -t $number -v 0.40
python util.py -d 0222 -m random -n 20000 -t $number -v 0.45
python util.py -d 0222 -m random -n 20000 -t $number -v 0.50
done
done 
exit 0
