void     next_line ();

void next_line(FILE *filep)

     /*
      *  Skips to next input line.
      */
{
  int dummy;
  while( (dummy=getc(filep)) != '\n');
}
