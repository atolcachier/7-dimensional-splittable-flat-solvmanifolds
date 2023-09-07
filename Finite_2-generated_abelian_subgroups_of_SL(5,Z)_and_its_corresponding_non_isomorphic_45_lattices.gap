#Script written in GAP to extract from the finite subgroups of GL(6,Z) the 2-generated abelian subgroups of SL(6,Z) which give rise to a lattice
#of a flat Lie group R^2 \ltimes \R^5 and distinguish the corresponding lattices
#
#Download cryst5 from https://drive.google.com/drive/folders/1E8PKEZ6mAM9oCUyJt-A7u76fJdZE3trH?usp=sharing and copy all the text from the .txt in GAP 
#(I extracted the .txt's from https://www.math.kyoto-u.ac.jp/~yamasaki/Algorithm/RatProbAlgTori/crystdat.html and gave them a format to copy into GAP without 
# changes)
#
2conj5:=[];
for i in [1..6079] do if Size(GeneratorsOfGroup(Group(cryst5[i])))=2 and IsAbelian(Group(cryst5[i]))=true and Determinant(GeneratorsOfGroup(Group(cryst5[i]))[1])=1
and Determinant(GeneratorsOfGroup(Group(cryst5[i]))[2])=1 then Append(2conj5, [cryst5[i]]);fi;od;
#
# Returns 105 subgroups. One can check that by computing
Size(2conj5);
#
# We construct the groups G_{A,B}:=Z^2 \ltimes_{A,B} \Z^4 and discard those with even rank of its abelianization #
F:=FreeGroup("a","b","c","d","e","t","s");
AssignGeneratorVariables(F);
C:=[];
D:=[];
for i in [1..105] do C[i]:=Matrix(2conj5[i][1])^-1;od;
for i in [1..105] do D[i]:=Matrix(2conj5[i][2])^-1;od;
G:=[];
for i in [1..105] do G[i]:=F/[Comm(a,b),Comm(a,c), Comm(a,d),Comm(a,e),Comm(b,c),Comm(b,d),Comm(b,e),Comm(c,d),Comm(c,e),Comm(d,e),Comm(t,s),
a^t*(a^C[i][1,1]*b^C[i][2,1]*c^C[i][3,1]*d^C[i][4,1]*e^C[i][5,1])^-1, 
b^t*(a^C[i][1,2]*b^C[i][2,2]*c^C[i][3,2]*d^C[i][4,2]*e^C[i][5,2])^-1, 
c^t*(a^C[i][1,3]*b^C[i][2,3]*c^C[i][3,3]*d^C[i][4,3]*e^C[i][5,3])^-1, 
d^t*(a^C[i][1,4]*b^C[i][2,4]*c^C[i][3,4]*d^C[i][4,4]*e^C[i][5,4])^-1, 
e^t*(a^C[i][1,5]*b^C[i][2,5]*c^C[i][3,5]*d^C[i][4,5]*e^C[i][5,5])^-1,
a^s*(a^D[i][1,1]*b^D[i][2,1]*c^D[i][3,1]*d^D[i][4,1]*e^D[i][5,1])^-1, 
b^s*(a^D[i][1,2]*b^D[i][2,2]*c^D[i][3,2]*d^D[i][4,2]*e^D[i][5,2])^-1, 
c^s*(a^D[i][1,3]*b^D[i][2,3]*c^D[i][3,3]*d^D[i][4,3]*e^D[i][5,3])^-1, 
d^s*(a^D[i][1,4]*b^D[i][2,4]*c^D[i][3,4]*d^D[i][4,4]*e^D[i][5,4])^-1, 
e^s*(a^D[i][1,5]*b^D[i][2,5]*c^D[i][3,5]*d^D[i][4,5]*e^D[i][5,5])^-1 ];od;
indexes:=[];
for i in [1..105] do if IsEvenInt(PositionNonZero(AbelianInvariants(G[i])))=true then Append(indexes, [i]); fi;od;
#
# Returns 45 subgroups, which are the following:. 
# 
lattices5:=[ [ [ [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ] ] ],
  [
      [ [ -1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 0, -1 ], [ 0, 0, 0, -1, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, -1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ -1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 0, -1 ], [ 0, 0, 0, -1, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, -1, 0, 0 ],
          [ 0, 0, 0, 0, 1 ], [ 0, 0, 0, 1, 0 ] ] ],
  [
      [ [ -1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 1, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, -1, 1 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, 0, -1, 1 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 1, 1, 0 ] ] ],
  [
      [ [ -1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, -1, 0, -1 ], [ 0, 0, -1, -1, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, 0, 1, -1 ],
          [ 0, 0, 1, 0, 1 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ],
          [ 0, -1, 0, 0, 0 ], [ 0, 0, -1, 0, 0 ] ],
      [ [ -1, 0, 0, 0, 0 ], [ 0, 0, 0, 1, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 1, 0, 0, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, -1 ], [ 0, 0, 0, -1, 0 ],
          [ 0, 0, -1, 0, 0 ], [ 0, -1, 0, 0, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1 ], [ 0, 0, 0, 1, 0 ],
          [ 0, 0, 1, 0, 0 ], [ 0, 1, 0, 0, 0 ] ] ],
  [
      [ [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ], [ 0, 0, 1, 0, 0 ],
          [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ] ],
      [ [ 0, 0, 0, 1, 0 ], [ 0, 0, -1, 0, 1 ], [ 0, 0, 1, 0, 0 ],
          [ 1, 0, 0, 0, 0 ], [ 0, 1, 1, 0, 0 ] ] ],
  [
      [ [ 0, -1, 0, 0, 0 ], [ -1, 0, 0, 0, 0 ], [ 0, 0, 0, -1, -1 ],
          [ 0, 0, -1, 0, -1 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ 0, 1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 1, 0, 1 ], [ 0, 0, -1, 1, 0 ] ] ],
  [
      [ [ 0, 1, 0, 0, 0 ], [ -1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 0, 0, -1 ],
          [ 0, 0, 1, 0, 1 ], [ 0, 0, 0, -1, 1 ] ],
      [ [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, 0, 1, -1 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, -1, -1, 0 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 0, 0, -1 ],
          [ 0, 0, 1, 0, 1 ], [ 0, 0, 0, -1, 1 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, -1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ] ] ],
  [
      [ [ 0, 1, 0, 0, 0 ], [ -1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, -1, 0, 0 ],
          [ 0, 0, 0, 0, 1 ], [ 0, 0, 0, 1, 0 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 1, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, -1, 0, 0, -1 ], [ 0, -1, 0, 0, 0 ] ],
      [ [ -1, 0, 0, 0, 0 ], [ 0, -1, -1, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 1, 0, -1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, -1, 0, 0 ], [ 0, 1, 0, 1, 1 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, -1, -1 ], [ 0, 0, -1, 1, 1 ],
          [ 0, 0, 0, 0, 1 ], [ 0, 0, 0, 1, 0 ] ] ],
  [
      [ [ 0, 0, 0, -1, 0 ], [ 1, 0, 0, 1, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, -1, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, 0, 0, 1 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 1, 0, 0 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 0, -1 ], [ 0, 0, 0, 1, 0 ] ],
      [ [ 0, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ 0, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 1, 0 ],
          [ 0, 0, -1, 0, -1 ], [ 0, 0, -1, 0, 0 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, -1, -1, 1 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 1, 1, 1, -1 ], [ 0, 1, 1, 1, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, -1 ], [ 0, 1, 0, 0, 1 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, -1, 0, 1 ] ] ],
  [
      [ [ 0, -1, -1, -1, 1 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 1, 1, 1, 1, -1 ], [ 1, 1, 1, 1, 0 ] ],
      [ [ 1, 1, 0, 0, 0 ], [ -1, 0, 0, 0, -1 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ -1, 0, 0, 0, 0 ] ] ],
  [
      [ [ 0, 0, -1, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 1, 0, 0, -1, 1 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ 1, -1, 0, -1, 1 ], [ 0, 1, 0, 1, -1 ], [ 0, 1, 1, 1, -1 ],
          [ 0, -1, 0, 0, 1 ], [ 0, 1, 0, 1, 0 ] ] ],
  [
      [ [ 0, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 0, -1 ], [ 0, 0, 0, 1, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ] ] ],
  [
      [ [ 0, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 1, 0 ],
          [ 0, 0, -1, 0, -1 ], [ 0, 0, -1, 0, 0 ] ],
      [ [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, -1, -1, -1, 0 ], [ 0, 1, 0, 0, 1 ],
          [ 0, 1, 1, 1, -1 ], [ 0, 1, 0, 1, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, -1, 0, 1 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 1, 1, 0, 0 ] ] ],
  [
      [ [ 0, 0, -1, -1, 1 ], [ -1, 0, 0, 0, -1 ], [ 0, 0, 1, 0, 0 ],
          [ 1, 1, 1, 1, -1 ], [ 0, 1, 1, 1, -1 ] ],
      [ [ 0, -1, -1, 0, 1 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 1, 1, 1, 0, 0 ] ] ],
  [
      [ [ 0, -1, -1, -1, 1 ], [ 0, 1, 0, 1, -1 ], [ 1, 1, 0, 0, 0 ],
          [ 0, -1, 0, 0, 1 ], [ 0, 1, 0, 1, 0 ] ],
      [ [ -1, 0, 0, 1, -1 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, -1, -1, 1 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, -1, 0, 0 ], [ 0, 1, 0, 0, 0 ],
          [ 0, 0, 0, 0, 1 ], [ 0, 0, 0, -1, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1 ], [ 0, 0, 0, 1, 0 ],
          [ 0, 0, 1, 0, 0 ], [ 0, 1, 0, 0, 0 ] ] ],
  [
      [ [ 0, 1, 0, 0, 0 ], [ -1, 0, 0, 0, -1 ], [ 0, 0, 0, 1, 0 ],
          [ 0, 0, -1, 0, 1 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ 0, 0, 1, 0, -1 ], [ 0, 0, 0, 1, -1 ], [ 1, 0, 0, 0, 1 ],
          [ 0, 1, 0, 0, 1 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ 1, 0, 0, 1, -1 ], [ 0, 0, 0, 1, 0 ], [ -1, 0, 0, -1, 0 ],
          [ 0, -1, 0, 0, 0 ], [ 1, 0, 1, 0, 0 ] ],
      [ [ 0, 0, -1, 1, 0 ], [ 0, 0, 0, 0, 1 ], [ 0, 0, 1, 0, 0 ],
          [ 1, 0, 1, 0, 0 ], [ 0, 1, 0, 0, 0 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ] ],
      [ [ 1, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ -1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 0, 1 ], [ 0, 0, 0, 1, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ], [ 0, -1, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 0, 1 ], [ 0, 0, 0, -1, 0 ] ],
      [ [ 1, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 1, 0 ],
          [ 0, 0, -1, -1, 1 ], [ 0, 0, 0, -1, 1 ] ],
      [ [ 1, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ] ],
      [ [ 1, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 1 ], [ 0, 0, 0, -1, 0 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 1 ], [ 0, 0, 0, 1, 0 ],
          [ 0, 0, 1, 0, 0 ], [ 0, 1, 0, 0, 0 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 1, -1 ],
          [ 0, -1, 1, 0, 0 ], [ 0, 0, 1, 0, 0 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 1 ], [ 0, 0, 0, -1, 0 ] ],
      [ [ -1, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ] ] ],
  [
      [ [ -1, 1, 0, 0, 0 ], [ -1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, 0, 1, 0 ],
          [ 0, 0, 0, 0, -1 ], [ 0, 0, -1, 0, 0 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, -1 ], [ 0, 0, 0, 1, 0 ] ],
      [ [ 1, -1, 0, 0, 0 ], [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, 0, 1, -1 ], [ 0, -1, 0, 0, 0 ],
          [ 0, 0, -1, 0, 1 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, -1 ], [ 0, 0, -1, 0, 1 ],
          [ 0, 1, -1, 0, 1 ], [ 0, 1, -1, 1, 1 ] ] ],
  [
      [ [ -1, 0, 0, 0, 0 ], [ 0, -1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 1 ], [ 0, 0, 0, -1, 0 ] ],
      [ [ 0, -1, 0, 0, 0 ], [ 1, -1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 0 ], [ 0, 0, 0, 0, -1 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, -1, 1 ], [ 0, 0, 0, -1, 0 ] ],
      [ [ 0, -1, 0, 0, 0 ], [ 1, -1, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ] ],
  [
      [ [ -1, 1, 0, 0, 0 ], [ -1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], [ 0, 0, 0, 0, -1 ],
          [ 0, 0, -1, 0, 0 ], [ 0, 0, 0, 1, 0 ] ] ],
  [
      [ [ 1, 0, 0, 0, 0 ], [ 0, 1, 0, 1, 0 ], [ 0, 1, 0, 1, 1 ],
          [ 0, -1, 1, -1, -1 ], [ 0, 0, 0, 1, 1 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, -1, -1, 0 ], [ 0, 0, 0, -1, -1 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 1, 0, 0, 0 ] ] ],
  [
      [ [ 0, 0, 0, 1, 0 ], [ 0, 1, 0, 1, 0 ], [ 0, 0, 1, -1, 0 ],
          [ -1, 0, 0, -1, 0 ], [ 0, 0, 0, 1, 1 ] ],
      [ [ 1, 0, 0, 0, 0 ], [ 0, 0, -1, 0, 0 ], [ 0, 0, 0, 0, -1 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 1, 0, 0, 0 ] ] ],
  [ [ [ 0, -1, -1, -1, 0 ], [ 0, 1, 1, 1, 1 ], [ -1, -1, 0, -1, 0 ],
          [ 0, 0, 0, 1, 0 ], [ 0, 0, -1, -1, 0 ] ],
      [ [ 0, 1, 1, 0, 1 ], [ 1, 0, -1, 0, -1 ], [ 0, 1, 1, 1, 1 ],
          [ 0, -1, 0, 0, -1 ], [ 0, 0, 0, 0, 1 ] ] ] ];;
#
# Next we have to check, for each of the 45 subgroups, if two generating sets give rise to isomorphic subgroups of the form $\Z^2\ltimes \Z^5$.
# With the following routine we check the order of the generating matrix as well as the order of the subgroup:
for i in [1..45] do
Print([Order(lattices5[i][1]),Order(lattices5[i][2]),Order(Group(lattices5[i][1],lattices5[i][2])),i], " "); od;

# The subgroups 1--9 correspond to groups of order 4 generated by both matrices of order 2.
# The subgroups 10--16 and 22--29 correspond to groups of order 8 generated by a matrix A of order 4 and a matrix B of order 2.
# The subgroups 17--21 correspond to groups of order 16 generated by both matrices of order 4.
# The subgroups 30--31 and 34--35 correspond to groups of order 12 generated by a matrix A of order 2 and a matrix B of order 6.
# The subgroups 32--33 correspond to groups of order 24 generated by a matrix A of order 4 and a matrix B of order 6.
# The subgroups 36--39 correspond to groups of order 18 generated by a matrix A of order 3 and a matrix B of order 6.
# The subgroup 40 corresponds to a group of order 36 generated by both matrices of order 6.
# The subgroups 41--45 correspond to groups of order 9 generated by both matrices of order 3.

# First we check which are the 2-generating sets of each subgroup:
          for i in [1..45] do
          for j in [0..Order(lattices5[i][1])-1] do
          for k in [0..Order(lattices5[i][1])-1] do
          for l in [0..Order(lattices5[i][2])-1] do
          for m in [0..Order(lattices5[i][2])-1] do
          if Order(Group(lattices5[i][1]^j*lattices5[i][2]^l,lattices5[i][1]^k*lattices5[i][2]^m))=Order(Group(lattices5[i])) then
          if j<k then Print([i,j,l,k,m], " "); 
          elif j=k then 
          if l<m then Print([i,j,l,k,m], " ");
          fi; fi; fi; od; od; od; od; od;

# We will use that the group \Z^2\ltimes_{A,B} \Z^5 is isomorphic to:
# \Z^2\ltimes_{B,A} \Z^5
#  \Z^2\ltimes_{A^{-1},B} \Z^5
#   \Z^2\ltimes_{A,AB} \Z^5.
# Basically, we can change {A,B} for {A^{-1},B} or {A,B^{-1}} and we can also perform the sequence of operations {A,B}->{A,A^k B}
# (or {A,B}->{B,AB^k}). We will try to list the powers [j,k,l,m] of the generating sets {A^j B^k, A^l B^m} in a way such that it is clear which operations
# we have to perform to get from one generator set to the following.
# For instance:
# For i in [1..9], in each case the outcome are the following powers: 
  [0,1,1,0], [0,1,1,1], [1,0,1,1],
# which corresponds to {A,B}->{B,A}, {A,B}->{B,AB} and {A,B}->{A,AB}.
#
# For i in [10..16], A^{-1}=A^3 and in each case the outcome is the following:

[ 0, 1, 1, 0 ] {A,B} -> {B,A},
[ 0, 1, 1, 1 ] {A,B} -> {B,AB}, 
[ 1, 0, 1, 1 ] {A,B} -> {A,AB}
         [ 1, 0, 2, 1 ] {A,AB} -> {A,A^2 B}, 
                   [ 1, 0, 3, 1 ] {A,A^2 B} -> {A,A^3 B}     
                   [ 2, 1, 3, 0 ] {A,A^2 B} -> {A^2 B, A^3}        (since A^3=A^{-1})
                   [ 2, 1, 3, 1 ] {A,A^2 B} -> {A^2 B, A^3 B}         

         [ 1, 1, 2, 1 ] {A,AB} -> {AB,A^2 B}, 
                   [ 1, 1, 3, 0 ] {AB,A^2 B} -> {AB,A^3}               

[ 0, 1, 3, 0 ] {A,B} -> {A^3,B} ->{B,A^3},   
                 [ 3, 0, 3, 1 ] {B,A^3} -> {A^3,A^3B}, 
                 [ 0, 1, 3, 1 ] {B,A^3} -> {B,A^3B},

# For i in [17..21], A^{-1}=A^3 and B^{-1}=B^3, and in each case the outcome is the following:

[ 0, 1, 1, 0 ] {A,B}->{B,A}
[ 0, 1, 1, 1 ] {A,B}->{B,AB}     
       [ 0, 1, 1, 2 ] {B,AB} -> {B,AB^2}  [ 0, 1, 1, 3 ] {B,AB^2} -> {B,AB^3} 
                                          [ 1, 2, 1, 3 ] {B,AB^2} -> {AB^2, AB^3} 
                                                    [ 1, 3, 2, 1 ]   {AB^2, AB^3} -> {A^2B, AB^3}
                                                                     [ 1, 3, 3, 0 ]  {A^2B, AB^3} -> {AB^3, A^3}
                                                                     [ 2, 1, 3, 0 ]  {A^2B, AB^3} -> {A^2B, A^3}
                                                    [ 1, 2, 2, 1 ]   {AB^2, AB^3} -> {AB^2, A^2B}
                                                                    [ 2, 1, 3, 3 ]   {AB^2, A^2B} -> {A^2B, A^3B^3}
                                                                    [ 1, 2, 3, 3 ]   {AB^2, A^2B} -> {AB^2, A^3B^3}
       [ 1, 1, 1, 2 ] {B,AB} -> {AB,AB^2} [ 1, 1, 2, 3 ] {AB,AB^2} -> {AB,A^2B^3}
                                                        [ 2, 3, 3, 0 ] {AB,A^2B^3} -> {A^2B^3, A^3}
                                          [ 1, 2, 2, 3 ] {AB,AB^2} -> {AB^2, A^2B^3}
                                                        [ 2, 3, 3, 1 ] {AB^2, A^2B^3} -> {A^2B^3, A^3B}
[ 1, 0, 1, 1 ] {A,B}->{A,AB}     
       [ 1, 0, 2, 1 ] {A,AB} -> {A,A^2B}   
                 [ 1, 0, 3, 1 ] {A,A^2B} -> {A,A^3B}
       [ 1, 1, 2, 1 ] {A,AB} -> {AB, A^2B} 
                  [ 1, 1, 3, 2 ] {AB,A^2B} -> {AB,A^3B^2}
                  [ 2, 1, 3, 2 ] {AB,A^2B} -> {A^2B,A^3B^2}
      [ 1, 1, 3, 0 ]  {A,AB} ->{AB,A^3}       
                
[ 0, 3, 1, 0 ] {A,B}->{B^3,A} 
        [ 0, 3, 1, 3] {B^3,A} ->{B^3, AB^3} 
                  [ 0, 3, 1, 2 ] {B^3,AB^3} ->{B^3,AB^2} 
                                [ 0, 3, 1, 1] {B^3,AB^2}->{B^3,AB}  
       [ 1, 0, 1, 3 ] {B^3, A} -> {A,AB^3} 
       [ 3, 0, 3, 3 ] {B^3, A} -> {A^3, A^3B^3}
                  [ 1, 0, 2, 3 ] {A,AB^3} -> {A,A^2B^3} 
                              [ 1, 0, 3, 3 ] {A,A^2B^3} -> {A,A^3B^3}
                              [ 2, 3, 3, 3 ] {A,A^2B^3} -> {A^2B^3, A^3B^3}
                  [ 1, 3, 2, 3 ] {A,AB^3} -> {AB^3, A^2B^3}
                               [ 2, 3, 3, 2 ] {AB^3, A^2B^3} -> {A^2B^3, A^3B^2}
                               [ 1, 3, 3, 2 ] {AB^3, A^2B^3} -> {AB^3, A^3B^2}
[ 0, 1, 3, 0 ] {A,B}->{B,A^3}  
       [ 0, 1, 3, 1 ] {B,A^3} -> {B,A^3B} 
                 [ 0, 1, 3, 2 ]  {B,A^3B} -> {B,A^3B^2} 
                             [ 0, 1, 3, 3 ]  {B,A^3B^2} -> {B,A^3B^3}
                             [ 3, 2, 3, 3 ]  {B,A^3B^2} -> {A^3B^2, A^3B^3}
                 [ 3, 1, 3, 2 ]  {B,A^3B} -> {A^3B,A^3B^2}
       [ 3, 0, 3, 1 ] {B,A^3} -> {A^3, A^3B} 
                  [ 2, 1, 3, 1 ] {A^3, A^3B} -> {A^2B, A^3B}
                                 [ 1, 2, 3, 1 ] {A^2B, A^3B} -> {AB^2,A^3B}

[ 0, 3, 3, 0 ] {A,B}->{B^3,A^3} 
       [ 0, 3, 3, 3 ] {B^3,A^3}->{B^3,A^3B^3} 
                  [ 0, 3, 3, 2 ] {B^3,A^3B^3}->{B^3,A^3B^2}  
                                [ 0, 3, 3, 1 ] {B^3, A^3B^2} - > {B^3,A^3B} 

# For i in [30..31] and i in [34..35] B^5=B^{-1} and in each case the outcome is the following:

[ 0, 1, 1, 0 ] [ 0, 1, 1, 1 ] [ 0, 1, 1, 2 ] [ 0, 1, 1, 3 ] [ 0, 1, 1, 4 ] [ 0, 1, 1, 5 ]  {B,A} -> {B, AB} -> {B, AB^2} -> .. -> {B,AB^5}
{B, A^k B} can be changed for {B^5, A^k B} so we obtain:
[ 0, 5, 1, 0 ] [ 0, 5, 1, 1 ] [ 0, 5, 1, 2 ] [ 0, 5, 1, 3 ] [ 0, 5, 1, 4 ] [ 0, 5, 1, 5 ]   

[ 1, 0, 1, 1 ] {A,B} -> {A,AB}
[ 1, 0, 1, 5 ] {A,B} -> {A,B^5} -> {A,AB^5}
[ 1, 1, 1, 2 ] {B,A} -> {B,AB} -> {AB,AB^2}
[ 1, 1, 1, 4 ] {AB,AB^2} ->{AB, AB^4}  (since AB^4=(AB^2)^-1)

[ 0, 3, 1, 1 ] {AB,AB^2} -> {B^3,AB}
[ 0, 3, 1, 4 ] {B^3,AB} -> {B^3,AB^4}
[ 0, 3, 1, 2 ] {B^3,AB^4} -> {B^3,AB^2}  (since AB^4=(AB^2)^-1)
[ 0, 3, 1, 5 ] {B^3,AB^2} -> {B^3, AB^5}

[ 1, 2, 1, 5 ] {B^3, AB^5} -> {AB^2, AB^5}
[ 1, 2, 1, 3 ] {AB^2, AB^5} -> {AB^2, AB^3}  (since (AB^2)^2*AB^5=AB^3)
[ 1, 4, 1, 5 ] {AB^2, AB^5} -> {AB^4, AB^5} (since AB^5=(AB^2)^-1)

[ 1, 3, 1, 4 ] {AB^2, AB^3} -> {AB^3, AB^4} (since AB^4=(AB^2)^-1)

# For i in [36..39] A^2=A^-1, B^5=B^-1 and in each case the outcome is the following:

[ 0, 1, 1, 0 ] [ 0, 1, 1, 1 ] [ 0, 1, 1, 2 ] [ 0, 1, 1, 3 ] [ 0, 1, 1, 4 ] [ 0, 1, 1, 5 ]  {B,A} -> {B,AB} -> {B,AB^2} -> ... -> {B,AB^5}
[ 0, 5, 1, 0 ] [ 0, 5, 1, 1 ] [ 0, 5, 1, 2 ] [ 0, 5, 1, 3 ] [ 0, 5, 1, 4 ] [ 0, 5, 1, 5 ]  (Changing B for B^5)
[ 0, 5, 2, 0 ] [ 0, 5, 2, 5 ] [ 0, 5, 2, 4 ] [ 0, 5, 2, 3 ] [ 0, 5, 2, 2 ] [ 0, 5, 2, 1 ] (Changing AB^k for its inverse)
[ 0, 1, 2, 0 ] [ 0, 1, 2, 5 ] [ 0, 1, 2, 4 ] [ 0, 1, 2, 3 ] [ 0, 1, 2, 2 ] [ 0, 1, 2, 1 ] (Changing again B^5 for B)

[ 1, 0, 1, 1 ] [ 1, 0, 2, 1 ] {A,B} -> {A,AB} -> {A, A^2B}
[ 2, 0, 1, 1 ] [ 2, 0, 2, 1 ] (Changing A for A^2)
[ 2, 0, 2, 5 ] [ 2, 0, 1, 5 ] (Changing A^kB for its inverse)
[ 1, 0, 2, 5 ] [ 1, 0, 1, 5 ] (Changing again A^2 for A) 

[ 0, 2, 1, 1 ] [ 0, 2, 1, 3 ] [ 0, 2, 1, 5 ]  {A,AB} -> {B^2, AB} (since A(AB)^2=B^2) -> {B^2, AB^3} -> {B^2, AB^5}
[ 0, 4, 1, 1 ] [ 0, 4, 1, 3 ] [ 0, 4, 1, 5 ]  (Changing B^2 for its inverse B^4) 
[ 0, 4, 2, 5 ] [ 0, 4, 2, 3 ] [ 0, 4, 2, 1 ]  (Changing AB^k for its inverse)
[ 0, 2, 2, 5 ] [ 0, 2, 2, 3 ] [ 0, 2, 2, 1 ]  (Changing again B^2 for B^4)

[ 1, 1, 1, 2 ] [ 1, 1, 2, 3 ] {B,AB} -> {AB,AB^2} -> {AB, A^2B^3}
[ 2, 5, 1, 2 ] [ 2, 5, 2, 3 ] (Changing AB for its inverse A^2B^5)
[ 2, 5, 2, 4 ] [ 2, 5, 1, 3 ] (Changing the second matrix for its inverse)
[ 1, 1, 2, 4 ] [ 1, 1, 1, 3 ] (Changing again A^2B^5 for AB)

 [ 1, 1, 2, 1 ] {A,AB} -> {AB,A^2B}
 [ 2, 5, 2, 1 ] (Changing AB for its inverse)
 [ 2, 5, 1, 5 ] (Changing A^2B for its inverse)
 [ 1, 1, 1, 5 ] (Changing again A^2B^5 for AB) 

 [ 1, 2, 2, 3 ] {AB,AB^2} -> {AB^2, A^2B^3}
 [ 2, 4, 2, 3 ] (Changing first matrix for its inverse)
 [ 2, 4, 1, 3 ] (Changing second matrix for its inverse)
 [ 1, 2, 1, 3 ] (Changing again first matrix for its inverse)

[ 1, 3, 1, 5 ] {B^2, AB^3} -> {AB^3, AB^5}
[ 2, 3, 1, 5 ] (Changing first matrix for its inverse) 
[ 2, 3, 2, 1 ] (Changing second matrix for its inverse)
[ 1, 3, 2, 1 ] (Changing again first matrix for its inverse)

[ 1, 4, 1, 3 ] {B^5,AB^4} -> {AB^4, AB^3} 
[ 2, 2, 1, 3 ] (Changing first matrix for its inverse)
[ 2, 2, 2, 3 ] (Changing second matrix for its inverse)
[ 1, 4, 2, 3 ] (Changing again first matrix for its inverse)

[ 1, 4, 1, 5 ] {B^5, AB^5} -> {AB^4, AB^5} 
[ 2, 2, 1, 5 ] (Changing first matrix for its inverse)
[ 2, 2, 2, 1 ] (Changing second matrix for its inverse)
[ 1, 4, 2, 1 ] (Changing again first matrix for its inverse)




# For i in [41..45], A^{-1}=A^2, B^{-1}=B^2, and in each case the outcome is the following:

[ 0, 1, 1, 0 ] {B,A}
[ 0, 1, 2, 0 ] {B,A} -> {B,A^2}
         [ 0, 1, 2, 1 ] {B,A^2} -> {B, A^2B}
                    [ 0, 1, 2, 2 ] {B,A^2B} -> {B, A^2B^2}

[ 0, 1, 1, 1 ] {B,A} -> {B,AB}
         [ 0, 1, 1, 2 ] {B,AB} -> {B,AB^2}
                   [ 1, 0, 1, 2 ] {B,AB^2} -> {A,AB^2}
                               [ 1, 0, 2, 2 ] {A,AB^2} -> {A,A^2B^2}
                                           [ 0, 2, 2, 2 ] {A,A^2B^2} -> {B^2, A^2B^2}
                                                         [ 0, 2, 2, 1 ] {B^2, A^2B^2} -> {B^2, A^2B} 
                                                                          [ 0, 2, 2, 0 ] {B^2, A^2B} -> {B^2, A^2}
                                                                                         [ 2, 0, 2, 2 ] {B^2, A^2} -> {A^2, A^2B^2}
                                                                                                       [ 1, 2, 2, 0 ] {A^2, A^2B^2} -> {AB^2, A^2}
                                                                                                                         [ 0, 2, 1, 2 ] {AB^2,A^2} -> {B^2, AB^2}
                                                                                                                                       [ 0, 2, 1, 1 ] {B^2, AB^2} -> {B^2, AB}
                                                                          [ 2, 0, 2, 1 ] {B^2, A^2B} -> {A^2, A^2B}
                                                                                         [ 1, 1, 2, 1 ] {A^2, A^2B} -> {AB, A^2B}
                                                                                         [ 1, 1, 2, 0 ] {A^2, A^2B} -> {AB, A^2}
                                                         [ 2, 1, 2, 2 ] {B^2, A^2B^2} -> {A^2B,A^2B^2}
                                           [ 0, 2, 1, 0 ] {A,A^2B^2} -> {B^2, A}
                                                                      
                               [ 1, 2, 2, 2 ] {A,AB^2} -> {AB^2, A^2B^2}
         [ 1, 1, 1, 2 ] {B,AB} -> {AB, AB^2}
[ 1, 0, 1, 1 ] {B,A} -> {A,AB}
         [ 1, 0, 2, 1 ] {A,AB} -> {A,A^2B}
 



for i in [36] do
          for j in [0..Order(lattices5[i][1])-1] do
          for k in [0..Order(lattices5[i][1])-1] do
          for l in [0..Order(lattices5[i][2])-1] do
          for m in [0..Order(lattices5[i][2])-1] do
          if Order(Group(lattices5[i][1]^j*lattices5[i][2]^l,lattices5[i][1]^k*lattices5[i][2]^m))=Order(Group(lattices5[i])) then
          if j<k then Print([j,l,k,m], " "); 
          elif j=k then 
          if l<m then Print([j,l,k,m], " ");
          fi; fi; fi; od; od; od; od; od;



# Next we distinguish the 45 lattices 
#
indexes2:=[];
for i in [1..Size(indexes)] do indexes2[i]:=Size(LowIndexSubgroupsFpGroup(G[indexes[i]],9));od;
for i in [1..Size(indexes2)] do
for j in [(i+1)..Size(indexes2)] do
if indexes2[i]=indexes2[j] then Print([i,j], " ");fi;od;od;
#
# Since there is nothing printed by this command, it means that 45 lattices are pairwise non-isomorphic.
#
# Here is the list of the 45 subgroups corresponding to these lattices:
