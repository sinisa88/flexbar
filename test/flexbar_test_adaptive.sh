
flexbar --source test_adaptive_right.fastq --target result_adaptive_right --format fastq-illumina1.5 --adapter-min-overlap 6 --adapter-seq AAAAAA  --min-readlength 10 --adapter-threshold 1 --adapter-trim-end RIGHT > /dev/null

a=`diff correct_result_adaptive_right.fastq result_adaptive_right.fastq`

#b=`diff correct_result_adaptive_right.fastq.omitted result_adaptive_right.fastq.omitted`

l1=`expr length "$a"`
#l2=`expr length "$b"`

if [ $l1 != 0 ]; then
echo "error testing mode adaptive-overlap, test 1"
echo $a
exit -1
else
echo "Test 1 OK"
fi


flexbar --source test_adaptive_left.fastq --target result_adaptive_left --format fastq-illumina1.5 --adapter-min-overlap 6 --adapter-seq AAAAAA  --min-readlength 10 --adapter-threshold 1 --adapter-trim-end LEFT > /dev/null

a=`diff correct_result_adaptive_left.fastq result_adaptive_left.fastq`

#b=`diff correct_result_adaptive_left.fastq.omitted result_adaptive_left.fastq.omitted`

l1=`expr length "$a"`
#l2=`expr length "$b"`

if [ $l1 != 0 ]; then
echo "error testing mode adaptive-overlap, test 2"
echo $a
exit -1
else
echo "Test 2 OK"
fi


echo ""

