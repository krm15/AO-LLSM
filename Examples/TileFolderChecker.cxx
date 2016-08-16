/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>

#include "itkDirectory.h"

int main ( int argc, char* argv[] )
{
  if ( argc < 3 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputFolder1 inputFolder2" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::Directory DirectoryType;

  std::string filename1, filename2;
  std::string searchString1 = "t.tif";
  std::string appendString = "_decon.tif";
  unsigned int numOfTIFs1(0), numOfTIFs2(0);

  DirectoryType::Pointer directory1 = DirectoryType::New();
  directory1->Load( argv[1] );
  unsigned int numOfFiles1 = directory1->GetNumberOfFiles();

  DirectoryType::Pointer directory2 = DirectoryType::New();
  directory2->Load( argv[2] );
  unsigned int numOfFiles2 = directory2->GetNumberOfFiles();

  for ( unsigned int m = 0; m < numOfFiles1; m++ )
  {
    filename1 = directory1->GetFile( m );
    if ( filename1.find( searchString1 ) != std::string::npos )
    {
      numOfTIFs1++;
    }
  }

  for ( unsigned int n = 0; n < numOfFiles2; n++ )
  {
    filename2 = directory2->GetFile( n );
    if ( filename2.find( appendString ) != std::string::npos )
    {
      numOfTIFs2++;
    }
  }

  std::cout << "Folder 1 has " << numOfTIFs1 << std::endl;
  std::cout << "Folder 2 has " << numOfTIFs2 << std::endl;

  for ( unsigned int m = 0; m < numOfFiles1; m++ )
  {
    //std::cout << "m: " << m << std::endl;
    filename1 = directory1->GetFile( m );
    if ( filename1.find( searchString1 ) != std::string::npos )
    {
      unsigned int lastindex = filename1.find_last_of(".");
      std::string filenameToBeFound = filename1.substr(0, lastindex) + appendString;

      //std::cout << "File to be found: " << filenameToBeFound << std::endl;

      bool found = false;
      for ( unsigned int n = 0; n < numOfFiles2; n++ )
      {
        filename2 = directory2->GetFile( n );
        if ( ( filenameToBeFound == filename2 ) && ( !found ) )
        {
          //std::cout << "File found: " << filename2 << std::endl;
          found = true;
          break;
        }
      }

      if ( !found )
      {
        std::cout << filenameToBeFound << std::endl;
      }
    }
  }

  return EXIT_SUCCESS;
}
