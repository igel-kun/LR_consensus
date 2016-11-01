
#pragma once

// a class storing for a fixed Mi and each vertex z of Si the size of a MAST between Mi and the subtree of Si rooted at z
// herein, i is in the range 0,..,p-1 and Mi is the subtree of the child of u_i that is not in the same centroid path as u_i
//NOTE: the sum of the sizes of all MiSi_Tables might be quadratic, but not all of it is going to be initialized
class MiSi_Table
{
  unsigned* mast_table;
  MyTree* const Mi;
  MyTree* const Si;

public:


  MiSi_Table(MyTree* const _Mi, MyTree* const _Si, const unsigned max_StId):
    mast_table(new unsigned[max_StId]), Mi(_Mi), Si(_Si)
  {}

  ~MiSi_Table()
  {
    delete Si;
    delete Mi;
    delete [] payload;
  }

  unsigned get_mast(const StId z) const
  {
  }

  unsigned compute_masts(const StId z)
  {
  }

  MyTree* get_Si() const
  {
    return Si;
  }
  

};


