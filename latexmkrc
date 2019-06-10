$pdflatex = 'pdflatex %O %S; pdfcrop --margins "20" %B.pdf %B-c.pdf; mv %B-c.pdf %B.pdf';
