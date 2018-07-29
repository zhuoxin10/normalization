BEGIN {
  dice = 0.0;
  md = 0.0;
  rmsd = 0.0;
  maxd = 0.0;
}
{ if (NR == 1) 
  {
    dice = $1;
  } else {
    md += $1;
    rmsd += $2;
    maxd += $3;
  }
}
END {
  printf "%02.6f %2.6f %2.6f %2.6f\n", dice, md/2.0, rmsd/2.0, maxd/2.0
}
