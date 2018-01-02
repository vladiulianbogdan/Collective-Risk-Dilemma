cpp="round2.out"

opt=1 #for the types of outputfiles 1=total ave, 2 round contrbutions per generation 12 or 21 for both...3 is special for old game!! dont use for sims it uses the same total files as new game

dumdum="dum.sh"

for initCont in -1 #0 1
do #-1 means rand initial state and any other number a homogenous pop so if you pick 0 then all population will contribute 0, if youi pick 2 they will contribute 2...mutations will change the initial state
for popsize0 in 100
do
for sel in 1 #selection intensity
do
for gamesize in 1000 #00 #how many games played in one generation
do
for mutype0 in 0.03 #0.001 0.01 0.1 0.03 #mutation rate
do
for gS0 in 1 #5 #10 #4 6 #10 20 #2 #4 #6 8 #group size
do
for rN in 1 4 #Round number
do
for randRound in -2 2 0 1 $rN # 0 -1 -2 3 4 #which round the risk will occur, 0 means all -1 non -2 one random, any Integer 1-rN means in that round
do
let endowMax=$rN\*\2;
for enom0 in 1 #endowment of each player at the begining of the game this is W in the paper
do
for enom1 in 4 #endowment of each player at the begining of the game this is W in the paper
do
for lostenom0 in 1 0.5 # this is alpha in the paper #part of the endowment lost with each risk round between 0-1
do
for lostenom1 in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 # this is alpha in the paper #part of the endowment lost with each risk round between 0-1
do
for maxPay0 in $enom0 #$enom0 #max amount one can pay per round
do
for maxPay1 in $enom1 #$enom0 #max amount one can pay per round
do
for probFnOPT0 in 3 #0=step curve 1=exp 1 2= fermi divided, 3 fermi function 4= exp 2 5 = piecewise
do
for prob0 in 1 #0 for 1 4 & 5, 1 for 0 2 & 3
do
for gamma0 in 10 #0.01 0.1 0.5 1 2 3 4 5 6 7 8 9 10 20 50 100 #lambda in the paper
do
for sigma1 in 0.15 #gaussian for threshold mutation
do
for sigma2 in -1 #gaussian for contributions mutation, -1 mean uniform mutation and other number implies the sigma for gausian
do
#let tD1=$enom0\+$enom2;
for target in 2 3 #this is 1 except for fermi or step curve opt 0 2 and 3
do 	

generatioMax=100000 #00#number of steps until the end of one simualtion

for num in 1 #{1..10} # iterations
	do
fname=$gamesize"-"$popsize0"-"$prob0"-"$sel"-"$gS0"-"$rN"-"$maxPay0"-"$maxPay1"-"$enom0"-"$enom1"-"$lostenom0"-"$lostenom1"-"$target"-"$probFnOPT0"-"$generatioMax"-"$gamma0"-"$sigma1"-"$sigma2"-"$mutype0"-"$mutype1"-"$randRound"-"$initCont
			subF="subM-"$fname

		
STR=$cpp" "$opt" 1 "$popsize0" "$sel" "$gamesize" "$generatioMax" "$mutype0" "$gS0" "$rN" "$maxPay0" "$maxPay1" "$enom0" "$enom1" "$lostenom0" "$lostenom1" "$target" "$prob0" "$probFnOPT0" "$gamma0" "$sigma1" "$sigma2" "$randRound" "$initCont" "$num #\${PBS_ARRAYID}

#echo $STR #call program with the parameters, must change this to mach your own termnial system
#./$STR
sname=$subF.sh
sed "s/dummy1/$STR/g;s/dum2/CR-round/g"  $dumdum > $sname
chmod a+x $sname
		
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done
done #
done
done
done
#done
#done
#done
#done
#done
#done


