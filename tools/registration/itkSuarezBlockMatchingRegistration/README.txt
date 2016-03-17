********************************************************************************

  This work has been made to be used against the libraries Insight Segmentation
  & Registration Toolkit 2.8.1.

  Program:    itkSuarezBlockMatchingRegistration
  Created by: Eduardo Suarez, Rafael Nebot
              Instituto Tecnologico de Canarias - Gobierno de Canarias
              http://www.itccanarias.org/
  Date:       2006/08/07


  License:

    Copyright (c) 2006 Instituto Tecnologico de Canarias - Gobierno de Canarias
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

     * Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.

     * Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

     * The name of the Instituto Tecnologico de Canarias - Gobierno de Canarias,
       may not be used to endorse or promote products derived from this software
       without specific prior written permission.
    
     * Modified source versions must be plainly marked as such, and must not be
       misrepresented as being the original software.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


  References:

  [1] Nonrigid Registration using Regularized Matching Weighted by
      Local Structure, E. Suarez, C.-F. Westin, E. Rovaris, J. Ruiz-Alzola
      MICCAI, Tokyo, Japan. 2002. Pages 581-589. 2002

  [2] Fast Entropy-Based Nonrigid Registration. E. Suarez, J. A. Santana,
      E. Rovaris, C.-F. Westin, J. Ruiz-Alzola. Computer Aided Systems Theory 
      (EUROCAST'03), Lecture Notes in Computer Science 2809.
      Las Palmas de Gran Canaria, Spain. Pages 607-615. 2003


  Acknowledgements:

    This work has been partially funded by project Torres Quevedo
    PTQ2004-1444 of the spanish goverment.

    The authors thanks to the National Alliance of Medical Image
    Computing (http://www.na-mic.org/) for his technical support in 
    the Insight Segmentation & Registration Toolkit.

    The authors also thank to Dan Blezek for his first approach to an
    implementation of this filter, to Luis Ibanez for his technical
    support.

********************************************************************************
