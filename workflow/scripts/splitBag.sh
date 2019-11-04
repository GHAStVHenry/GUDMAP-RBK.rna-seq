#!/bin

for i in $(ls -d Replicate_*)
do
rsync -r $1/ ${i} --exclude=fetch.txt
zip -r ${i}.zip ${i}
done