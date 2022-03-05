! This program is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
! 
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
!----------------------------------------------------------------------------
! Created By  : Daniel O. Perez
! Created Date: between 2010 and 2015 
! email: perez.daniel.omar@gmail.com
! ---------------------------------------------------------------------------


module count_file_lines_mod

  use precision, wp=>dp

contains
  
  subroutine count_file_lines(file_in,lines_out)

    !esta subrutina cuenta las lineas de un archivo
    !si el archivo no fue abierto, se abre y se cierra, sino se deja abierto
    
    implicit none
    character(len=*), intent(in)::file_in
    integer, intent(out)::lines_out
    integer::myunit,io
    real(kind=wp)aux

    !veo si el archivo está abierto
    inquire(file=file_in, number=myunit)

    if(myunit.le.0)then !si no está abierto
       lines_out=0
       open(newunit=myunit,file=file_in,action='read')
       do
          read(myunit,*,iostat=io)aux
          if(io/=0)exit
          lines_out=lines_out+1
       end do
       close(myunit)
    else !si está abierto
       lines_out=0
       do
          read(myunit,*,iostat=io)aux
          if(io/=0)exit
          lines_out=lines_out+1
       end do
       rewind(myunit)
    end if
       
    
  end subroutine count_file_lines

  
end module count_file_lines_mod
