SHELL=/bin/sh

# default
PCx:: 
	cd SRC;       ${MAKE}

PCx_NgPeyton:: 
	cd SRC;       ${MAKE} PCx_NgPeyton

PCx_wssmp:: 
	cd SRC;       ${MAKE} PCx_wssmp

mex::
	cd mex; 	${MAKE} mex

PCx_mysolver:: 
	cd SRC;       ${MAKE} PCx_mysolver

clean:
	/bin/rm -f PCx *.log ./Ng-Peyton/*.[ao] ./SRC/*.o 
	/bin/rm -f SRC/libPCx.a mex/*.mex* mex/*.[ao] F2C/*.[ao]
