ori_svg=$1
crop_svg=$2

cat $ori_svg | grep "path style" | sed 's/^.*M \([0-9]*\.\?[0-9]* [0-9]*\.\?[0-9]*\) .*$/\1/g' | sed '1d' > ./original.coord
cat $crop_svg | grep -P "^\t<path" | sed 's/^.*M\([0-9]*\.\?[0-9]*,[0-9]*\.\?[0-9]*\).*$/\1/g' | sed 's/,/ /g' | sed '1d' > ./cropped.coord

source activate base
conda activate daily
python /media/biogenger/D/scripts/GZLAB_ST_PIPELINE/filter_pixel.py \
--coorR=./original.coord \
--coorC=./cropped.coord \
--output=./manual.remove.txt
