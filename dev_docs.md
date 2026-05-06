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
a vector of `IBD_segment`s called `IBDs`, one for each sample but not the archaic.
It does this by using the SampleMapper class which handles user-specified sample subsets.

Initializing an `IBD_segment` requires a "name", "threshold", "pool", and
"exclusive" flag. These are all simple except for "pool" which I'm still
trying to figure out.
- "name" is the sample name
- "threshold" is the LOD-threshold specified by args
- "exclusive" is a bool for using `(start, end)` or `(start, end]`
- "pool"  I'm still figuring out but it gets passed by reference
  to all the `IBD_segments` so I think it's a single shared obj

I'm surprised that `IBDs.emplace_back(sample, threshold, &pool, exclusive_end);`
in `IBD_Collection::initialize` is valid since it seems to be implicitly creating
`IBD_segment` objects and adding them to the IBDs vector in a single command. Didn't
know this was possible with `c++`.

The rest of main.cc proceeds by alternating calls to GenotypeReader::update()
and `IBD_collection::update()`.

`GenotypeReader::update()` reads the next line of the genotype file,
checks if the locus is within the mask file then calculates does more
processing with `GenotypeReader::process_line_buffer()`. This includes
finding the frequency at this locus using `Genotype_Reader::find_frequency()`
which uses only modern genotypes, not the archaic. Then it uses a `LODCalculator`,
which is a member of the `GenotypeReader` class, to calculate the lodscore of
this locus for each individual. Since there are only a handful of possible lodscores
per locus because given the allele frequency, it only depends on the modern/ancient genotypes.
That's why there is a `calculator.update_lod_cache` and then a `calculator.calculate_lod` because
IBDmix is caching the scores to avoid recalculating the same thing over and over.
These results are stored in the `lod_scores` member.

`IBD_collection::update()` takes as input the `GenotypeReader` and calls
`IBD_segment::add_lod()` on each of the `IBD_segments` in the `IBDs` vector.
This gets its info from the `GenotypeReader` `lod_scores` which were calculated
in the `GenotypeReader::update()` function.

This call to `IBD_segment::add_lod()` calls `IBD_sgement::add_node` which uses
the `IBD_pool->get_node` function to create(?) a new `IBD_node` and "push"
it onto the `IBD_stack` which is named `segment`. I think the `IBD_pool` is
being smart about memory and is re-using prior `IBD_Node` memory allocations, but I need
to look more into this. NOTE that if this is true, `IBD_pool` might make it complicated
to parallelize ibdmix, or `IBD_pool->get_node` would need to have atomic operations? It
looks like each `IBD_Segment` has it's own `IBD_pool`.

Any, this `IBD_Node` get's pushed onto `IBD_Segment`'s `segment`
member, which is an `IBD_Stack`.  The `IBD_stack` is a singly-linked list of
`IBD_Node` which is a data struct:
```
struct IBD_Node {
  double cumulative_lod, lod;
  uint64_t position;
  unsigned char bitmask;
  IBD_Node *next;
};
```
Note that `IBD_Node` does not need a `chrom` field because all nodes
in the same stack are guaranteed to be on the same `chrom`.

During the `IBD_segment::add_node`, a new node of the current position
becomes the new head of the singly-linked list of `IBD_Node`s that make
up the `IBD_stack`. The `cumulative_lod` of the new node is the sum of
it's own lod and the lod of the prior top.

After the `IBD_node` is pushed to `IBD_stack`, there are some logically
cases to consider such as:
* We only have one node on the stack, so we're starting a new segment
  -> We have to re-initialize the recorders to set all their values to 0
  -> Then record this new locus for all recorders

* `segment.topIsNewMax()` the best cumulative lod score is now the just-added node.
  -> NOTE I'm not very sure about these descriptions
  -> We set the end of the segment to be the new top which means that
     we've expanded the length of the segment because this is a "good" node
     that we just added.
  -> Then we run `pool->reclaim_segment(&segment)` but I'm not sure why.
     This calls `IBD_Stack::getSegmentFrom(IBD_Stack *other)` which I'm
     getting confused about. Seems like it's joining the pool stack and the segment stack?

* `segment.reachedMax()` means that by adding this last node, the cumulative lod
   at the top is now negative.
   -> We write-out the best segment from start to end
   -> Then do something with unprocessing the segment which involves
      reversing the linked list
   -> Then we `pool->reclaim_stack(&segment)` which is somehow returning nodes
      back to the pool?
   -> Finally it looks like we are cleaning up the unprocessed segment?

Side note that the `pool` member of `IBD_Pool` is actually an `IBD_Stack` itself,
which I guess makes sense since I think `IBD_Pool` is trying to keep track of
already allocated nodes for reuse?

Main take-aways:
> Main loop alternates between (a) fetching the next line of the genotype file,
  calculating the lod per sample and (b) adding these lods to IBD_collection.
> A single `IBD_collection` is created which has a vector of `IBD_segments`,
> one for each sample. These `IBD_segments` each have an `IBD_stack` which is
> a singly-linked list of `IBD_node`s.

Things I still don't understand
- The purpose of the `IBD_pool`
- How the `IBD_stack` is creating the introgressed regions 