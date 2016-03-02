export LD_LIBRARY_PATH=../../RLCS/lib:$LD_LIBRARY_PATH

#wc -l MSS_Bhairavi2_MSSConcertIX.wav.spl.mel.441f0 | cut -d' ' -f1 > temp
#cut -d' ' -f2 MSS_Bhairavi2_MSSConcertIX.wav.spl.mel.441f0 >>  temp

../bin/get_regions_of_interest ../examples/temp2 100 50 ../examples/

