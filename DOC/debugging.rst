Debugging
=========

For debugging PCx you can use GDB [GDB]_.

A very simple example of the use of GDB can be found below. ::

    $ gdb PCx
    GNU gdb (GDB) 7.5.1
    Copyright (C) 2012 Free Software Foundation, Inc.
    License GPLv3+: GNU GPL version 3 or later
    <http://gnu.org/licenses/gpl.html>
    This is free software: you are free to change and redistribute it.
    There is NO WARRANTY, to the extent permitted by law.  Type "show copying"
    and "show warranty" for details.
    This GDB was configured as "x86_64-unknown-linux-gnu".
    For bug reporting instructions, please see:
    <http://www.gnu.org/software/gdb/bugs/>...
    Reading symbols from $HOME/PCx/PCx...done.
    (gdb) break main.c:70
    Breakpoint 1 at 0x401386: file main.c, line 70.
    (gdb) run
    Starting program: $HOME/PCx/PCx
    warning: Could not load shared library symbols for linux-vdso.so.1.
    Do you need "set solib-search-path" or "set sysroot"?

    Breakpoint 1, main (argc=1, argv=0x7fffffffe818) at main.c:70
    70         printf("\n******** PCx version 1.1 (Nov 1997) ************\n\n");
    (gdb) c
    Continuing.

    ******** PCx version 1.1 (Nov 1997) ************

    Usage:
            $HOME/PCx/PCx mpsfile
    [Inferior 1 (process 4757) exited with code 01]
    (gdb) quit

Tips
----

For print a structure you can use ::

    (gdb) p *strucuture

For print a array you can use ::

    (gdb) p *array@len

.. rubric:: Refereces

.. [GDB] Free Software Foundation. GDB: The GNU Project Debugger. http://www.gnu.org/software/gdb/.
