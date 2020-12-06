URL="http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp"

cd data/

while read line; do
file=$(echo $line | tr -d '\r')
  wget $URL/$file
done < instance_list.txt

cd -