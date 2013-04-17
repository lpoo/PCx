file ../../PCx
b Ng-Peyton.c:282
set loggin overwrite on
run maros0
def psfinit
set loggin on
p dimension
p nonzeros
p *Factor.Perm@dimension
p *Factor.InvPerm@dimension
p *ColumnCount
#p *Factor.NonzerosL
p NgPeyton.NumCompressedCols
set loggin off
end
set loggin file maros0_sfinit_before.txt
psfinit
n
set loggin file maros0_sfinit_after.txt
psfinit
run afiro
set loggin file afiro_sfinit_before.txt
psfinit
n
set loggin file afiro_sfinit_after.txt
psfinit
q
