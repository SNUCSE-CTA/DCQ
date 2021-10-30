num=4000
for d in 1 2
do
	./DCQ -m all -f data/hprd_D${num}_d${d}.igraph -q query/hprd_d${d}.igraph
done

d=2
for num in 1000 2000 4000 8000
do
	./DCQ -m all -f data/hprd_D${num}_d${d}.igraph -q query/hprd_d${d}.igraph
done
