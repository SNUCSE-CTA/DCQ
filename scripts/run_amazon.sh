num=4000
for d in 1 2 3
do
	./DCQ -m all -f data/amazon_D${num}_d${d}.igraph -q query/amazon_d${d}.igraph
done

d=2
for num in 1000 2000 4000 8000
do
	./DCQ -m all -f data/amazon_D${num}_d${d}.igraph -q query/amazon_d${d}.igraph
done
