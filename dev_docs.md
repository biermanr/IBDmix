March 17th 2026
---------------

I'm walking through IBDmix source code carefully now that I have
a slightly better understanding of (a) C++ (b) IBDmix itself and
(c) a desire to see if IBDmix can be parallelized (d) a desire
to see if IBDmix can use bcftools for random access.

I think I've done this at least once before, but as I look through it
now it all seems new, so I guess it didn't stick. I'll write my notes here.


So main.cc is what's used to create the `ibdmix` executable which is what
I'm most interested in. It starts by creating a `GenotypeReader` and then
creating an `IBD_collection`.

`IBD_collection::initialize` takes as input the `GenotypeReader` and creates
a vector of `IBD_segment`s called `IBDs`, one for each sample, including the archaic.

Initializing an `IBD_segment` requires a "name", "threshold", "pool", and
"exclusive" flag. These are all simple except for "pool" which I'm still
trying to figure out.
- "name" is the sample name
- "threshold" is the LOD-threshold specified by args
- "exclusive" is a bool for using (start, end) or (start, end]
- "pool"  I'm still figuring out but it gets passed by reference
  to all the `IBD_segments` so I think it's a single shared obj

Next during the `IBD_segment::initialize()` it sets up all the additional
"Recorders" such as the one used by `more-stats`.

The rest of main.cc proceeds by alternating calls to GenotypeReader::update()
and `IBD_collection::update()`.

`GenotypeReader::update()` reads the next line of the genotype file,
calculates the LOD score (of each sample?) and updates it's values.

`IBD_collection::update()` takes as input the `GenotypeReader` and calls
`IBD_segment::add_lod()` on each of the `IBD_segments` in the `IBDs` vector.

This call to `IBD_segment::add_lod()` calls `IBD_sgement::add_node` which uses
the `IBD_pool->get_node` function to create(?) a new `IBD_node` and "push"
it onto the `IBD_stack` which is named `segment`.

The `IBD_stack` is a singly-linked list of `IBD_Node` which is a data struct:
```
struct IBD_Node {
  double cumulative_lod, lod;
  uint64_t position;
  unsigned char bitmask;
  IBD_Node *next;
};
```

During the `IBD_segment::add_node`, a new node of the current position
becomes the new head of the singly-linked list of `IBD_Node`s that make
up the `IBD_stack`. There are important conditionals that I still need
to understand such as `topIsNewMax` and `reachedMax` which determine
what should happen to the `IBD_stack` now that the new node has been added.
It looks like there is logic for writing out a stack and clearning it when
we've come to the end of an introgressed segment presumably.

Ok, that's all for now. I have a better understanding of the "class hierarchy"
which goes:


> A single `IBD_collection` is created which has a vector of `IBD_segements`,
> one for each sample. These `IBD_segments` each have an `IBD_stack` which is
> a singly-linked list of `IBD_node`s.

Things I still don't understand
- How the `IBD_pool` comes into play
- How the `IBD_stack` is creating the introgressed regions 
- Other things I'm sure







