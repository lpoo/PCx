file ../../PCx
b Ng-Peyton.c:256
set loggin overwrite on
def pordmmd
set loggin on
p dimension
p *NgPeyton.pSuperNodeCols@dimension
p *NgPeyton.SuperNodeRows@NgPeyton.pSuperNodeCols[dimension]-1
p *Factor.InvPerm@dimension
p *Factor.Perm@dimension
p NgPeyton.NumCompressedCols
set loggin off
end
run maros0
set loggin file maros0_ordmmd_before.txt
pordmmd
n
set loggin file maros0_ordmmd_after.txt
pordmmd
run afiro
set loggin file afiro_ordmmd_before.txt
pordmmd
n
set loggin file afiro_ordmmd_after.txt
pordmmd
q
