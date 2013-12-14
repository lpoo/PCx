cc      program transforma
c      DOUBLE PRECISION A(6),b(3)
c      INTEGER*4 col(4),row(6),n,nnulos,iter
c      iter=1
c      n=3
c      nnulos=6
c      A(1)=3.
c      A(2)=1.
c      A(3)=-1.
c      A(4)=3.
c      A(5)=1.
c      A(6)=2.
c      col(1)=1
c      col(2)=4
c      col(3)=6
c      col(4)=7
c      row(1)=1
c      row(2)=2
c      row(3)=3
c      row(4)=2
c      row(5)=3
c      row(6)=3
c      b(1)=1.
c      b(2)=2.
c      b(3)=3.
c      call converte(iter,n,nnulos,row,col,A,b)
c      stop
c      end

C**********************************************************************
C     ***  Este programa recebe como entrada um  sitema Ax=b gerado pelo
C     ***  codigo PCx e armazena-o no formato Harwell-Boeing
C**********************************************************************


      SUBROUTINE transforma(iter,n,nnulos,row,col,A,b,solution)

c      SUBROUTINE TRANSFORMA(iter,n,nnulos,row,col,A,b,solution)
      DOUBLE PRECISION A(*),b(*),solution(*)
      INTEGER*4 col(*),row(*),n,nnulos,iter

C      parametros de entrada:
C      iter: iteracao do PCx
C      n :dimensao da matriz A
C      nnulos: no de elementos nao nulos de A
C      row: ponteiro para linhas de A
C      col: ponteiro para colunas de A
C      A: valores dos coeficientes armazenado por colunas
C      b: vetor dos termos independentes


C     parametros de saida
C     arquivo com o sistema no formato Harwell-Boeing


C variaveis locais
      integer
     X        TOTCRD,PTRCRD,INDCRD,VALCRD,RHSCRD,
     X        NROW,NCOL,NNZERO,NELTVL,NRHS,NRHSIX,
     X        aux,nlx,uss,i1,i2



      CHARACTER TITLE*72, KEY*8, MXTYPE*3,RHSTYP*3,conca*30,
     X          PTRFMT*16,INDFMT*16,VALFMT*20,RHSFMT*20



      iter=iter+1
      uss=2

c     transforma numero em caractere. Utiliza codigo ascii
      i2 = mod(iter,10)
      i1 = (iter-i2)/10
      print*,i1, ' ',i2
      conca = 'qap15_'//char(48+i1) // char(48+i2) // '.rsa'
      print *, ' conca', ' ', conca



      open(unit=uss,file=conca,status='new')

c     conta numero de linhas que serao necessarias p/ armazenar
c     os dados do arquivo
c     matriz A
      aux=nnulos-int(nnulos/4)*4
      if (aux.eq.0) then
         VALCRD=int(nnulos/4)
      else
         VALCRD=int(nnulos/4)+1
      endif

c     ponteiros p/ linhas de A
      aux=nnulos-int(nnulos/10)*10
      if (aux.eq.0) then
         INDCRD=int(nnulos/10)
      else
         INDCRD=int(nnulos/10)+1
      endif

c     ponteiros p/ colunas
      aux=n+1-int((n+1)/10)*10
      if (aux.eq.0) then
         PTRCRD=int((n+1)/10)
      else
         PTRCRD=int((n+1)/10)+1
      endif

c     solx,bx
      aux=n-int(n/4)*4
      if (aux.eq.0) then
         nlx=int(n/4)
      else
         nlx=int(n/4)+1
      endif


c     numero total de linhas
      TOTCRD=PTRCRD+INDCRD+VALCRD+2*NLX


      RHSCRD=nlx+nlx
      MXTYPE='RSA'
      NROW=n
      NCOL=n
      NNZERO=nnulos
      NELTVL=0
      PTRFMT='(10I8)'
      INDFMT='(10I8)'
      VALFMT='(4E18.8)'
      RHSFMT='(4E18.8)'
      RHSTYP='FBX'
      NRHS=1
      NRHSIX=n

      TITLE='Problema1'


c     escreva no arquivo HEADER BLOCK
      WRITE(uss,100),TITLE,ITER,
     X               TOTCRD,PTRCRD,INDCRD,VALCRD,RHSCRD,
     X               MXTYPE,NROW,NCOL,NNZERO,NELTVL,
     X               PTRFMT,INDFMT,VALFMT,RHSFMT
 100  format(A72,I8 / 5I14 / A3,11X, 4I14/ 2A16,2A20)
      WRITE(uss,101),RHSTYP,NRHS,NRHSIX
 101  format(A3,11X,2I14)
c     escreve no arquivo vetor de ponteiros p/ colunas de ACOEFF
      WRITE(uss,PTRFMT),(col(i),i=1,NCOL+1)
c     escreve no arquivo vetor de linhas p/ ACOEFF
      WRITE(uss,INDFMT),(row(i),i=1,NNZERO)
c     escreve no arquivo coeficientes nao nulos de ACOEFF
      WRITE(uss,VALFMT),(A(I),I=1,NNZERO)
c     escreve no arquivo os coeficientes b
      WRITE(uss,RHSFMT),(b(i),i=1,NCOL)
c     escreve no arquivo solucao exata
      WRITE(uss,RHSFMT),(solution(i),i=1,NCOL)

      close(unit=uss,iostat=iost,status='keep')


      RETURN
      END
