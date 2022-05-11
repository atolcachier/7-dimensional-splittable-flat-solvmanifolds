#Script written in GAP to extract from the finite subgroups of GL(6,Z) the 2-generated abelian subgroups of SL(6,Z) which give rise to a lattice
of a flat Lie group R^2 \ltimes \R^5 and distinguish the corresponding lattices /*

/*Download cryst5 from https://drive.google.com/drive/folders/1E8PKEZ6mAM9oCUyJt-A7u76fJdZE3trH?usp=sharing and copy all the text from the .txt in GAP 
(I extracted the .txt's from https://www.math.kyoto-u.ac.jp/~yamasaki/Algorithm/RatProbAlgTori/crystdat.html and gave them a format to copy into GAP without changes) /*

2conj5:=[];
for i in [1..6079] do if Size(GeneratorsOfGroup(Group(cryst5[i])))=2 and IsAbelian(Group(cryst5[i]))=true and Determinant(GeneratorsOfGroup(Group(cryst5[i]))[1])=1
and Determinant(GeneratorsOfGroup(Group(cryst5[i]))[2])=1 then Append(2conj5, [cryst5[i]]);fi;od;

/* Returns 105 subgroups. One can check that by computing
Size(2conj5);

/* We construct the groups G_{A,B}:=Z^2 \ltimes_{A,B} \Z^4 and discard those with even rank of its abelianization */
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

/* Returns 45 lattices. Next we distinguish these lattices */

indexes2:=[];
for i in [1..Size(indexes)] do indexes2[i]:=Size(LowIndexSubgroupsFpGroup(G[indexes[i]],9));od;
for i in [1..Size(indexes2)] do
for j in [(i+1)..Size(indexes2)] do
if indexes2[i]=indexes2[j] then Print([i,j], " ");fi;od;od;

/* Since there is nothing printed by this command, it means that 45 lattices are pairwise non-isomorphic.*/

