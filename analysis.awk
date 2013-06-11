# Only work for two specifications file.
BEGIN{
    TF=0;  # Total number of fail
    LF=0;  # Left number of fail
    RF=0;  # Right number of fail
    ZLT=0;  # Nonzeros <
    ZEQ=0;  # Nonzeros ==
    ZGT=0;  # Nonzeros >
    ILT=0;  # Iters <
    IEQ=0;  # Iters ==
    IGT=0;  # Iters >
    TLT=0;  # Time <
    TEQ=0;  # Time ==
    TGT=0;  # Time >
}
{
    if ($2 != 0 || $6 != 0){
        TF = TF + 1;
        if ($2 != 0) LF = LF + 1;
        if ($6 != 0) RF = RF + 1;
    }
    else{
        if ($3 < $7) ZLT = ZLT + 1;
        if ($3 == $7) ZEQ = ZEQ + 1;
        if ($3 > $7) ZGT = ZGT + 1;
        if ($4 < $8) ILT = ILT + 1;
        if ($4 == $8) IEQ = IEQ + 1;
        if ($4 > $8) IGT = IGT + 1;
        if ($5 < $9) TLT = TLT + 1;
        if ($5 == $9) TEQ = TEQ + 1;
        if ($5 > $9) TGT = TGT + 1;
    }
}
END{
    print "Total number of fail: " TF;
    print "Left number of fail: "  LF;
    print "Right number of fail: " RF;
    print "Nonzeros <: "           ZLT;
    print "Nonzeros ==: "          ZEQ;
    print "Nonzeros >: "           ZGT;
    print "Iters <: "              ILT;
    print "Iters ==: "             IEQ;
    print "Iters >: "              IGT;
    print "Time <: "               TLT;
    print "Time ==: "              TEQ;
    print "Time >: "               TGT;
}
