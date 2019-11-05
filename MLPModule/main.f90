Program nGAP
integer              :: imode, j
character(100)       :: line
integer              :: ios
imode = 1
open(2233, file = 'control')
DO while (.true.)
    read(2233,'(A100)', iostat=ios) line
    if (ios .ne. 0) exit
    if (index(line, 'lrelax') .ne. 0) then
        j = index(line, '=')
        read(line(j + 1:),*) imode
    endif
ENDDO
print*, 'imode',imode

call mlp_main(imode)
END PROGRAM
