        integer, parameter :: n = 369007
        integer nodes, fid, a
        integer, dimension(:) , allocatable :: roi
        CHARACTER*1024 filename
        filename = "baseline_ASI2_m1_nodes.txt"
        nodes = 0
        fid = 1 
        OPEN (fid, FILE=TRIM(fileName))
        do while(.true.)
           
        read(fid,*,END=111) 
           nodes = nodes + 1 
        end do
 111    rewind(fid)
        print * , nodes
        close(fid)
        allocate (roi(nodes))
        fid = 1
        OPEN (fid, FILE=TRIM(fileName))
        do k = 1, nodes
        
        read(fid,*) roi(k)
        
        end do

        close(fid)
    
    
        open(1,file="tagFile_bl_asi6_p3")
        write(1,*) n, 1
        
        do i=1, n
            
            if( ANY(roi .EQ. i)) then
                         
               write(1,*) 1
            else
               write(1,*) 0
            end if
        end do
        close(1)
        end
