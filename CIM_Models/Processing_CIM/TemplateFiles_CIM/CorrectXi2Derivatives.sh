#! /bin/tcsh

# This file is used to correct the Xi2 derivatives for the apex to remain circular 


foreach NODE (1) 
  awk -f CorrectXi2DerivativesSub.awk -v node1=$NODE temp.ipnode > test.ipnode
  echo "Updated node $NODE"
end

foreach NODE (18) 
  awk -f CorrectXi2DerivativesSub.awk -v node1=$NODE test.ipnode > testV2.ipnode
  echo "Updated node $NODE"
end

mv testV2.ipnode temp.ipnode
rm test.ipnode
