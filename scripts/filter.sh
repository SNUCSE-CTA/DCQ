G=$1
N=$2
D=$3

DATA_DIR="./data"
QUERY_DIR="./query"

output=$(./DCQ -m filter -f ${DATA_DIR}/${G}_D${N}_d${D}.igraph -q ${QUERY_DIR}/${G}_d${D}.igraph 2>&1)
x=$?

echo "${output}"

exit ${x}

