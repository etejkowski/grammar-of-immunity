# What this study asked, and what we found

*A plain-language summary. No background in immunology or computing assumed.*

Erick Tejkowski and Maria Elisa Paredes · Fairview Heights, Illinois

---

## The problem

Your immune system carries billions of T cells. Each one has a receptor on its
surface that recognizes one specific fragment of an invader — a virus protein, a
bacterial protein, a piece of a tumor. Think of the receptor as a key and the
fragment as a lock. When a key fits a lock, that T cell multiplies and attacks.

If we could look at a receptor's sequence and predict which fragment it opens, a
great deal would follow: vaccines designed with more confidence, cancer therapies
matched to a patient's own T cells, blood tests that read someone's infection
history from their immune repertoire.

We cannot do this yet. Computer models predict reasonably well for fragments they
were trained on, and collapse to roughly coin-flip accuracy for fragments they
have never seen. That second case is the one that matters, because a new virus is
by definition a fragment nobody trained on.

## The idea we tested

The receptors are not designed; they are assembled. Your body builds each one by
picking two inherited gene pieces — a beginning and an end — and stitching them
together with a short run of essentially random letters in the middle. The
critical, lock-touching part of the receptor is therefore made of three parts: an
inherited front, a random middle, and an inherited back.

In 1984, immunologist Niels Jerne suggested in his Nobel lecture that the immune
system works something like a language, with its own grammar. In 2024, a research
group made that concrete: they argued that current computer models fail partly
because they chop receptor sequences into meaningless fragments, the way you
might cut the word *unhappiness* into *unha*, *ppin*, *ess*. A linguist would
cut it into *un-*, *happy*, *-ness* — the parts that actually carry meaning. They
proposed doing the same for receptors, and invited someone to try it.

We tried it. Cut each receptor at its real seams — inherited front, random middle,
inherited back — and see whether models built from those meaningful parts predict
better, especially for fragments they have never seen.

## What we found

**The structure is real.** The "random" middle is not fully random. We could tell
a genuine middle from the same letters shuffled into a different order, which
means the order follows rules. Knowing which inherited pieces a receptor used
made that easier still.

**But it doesn't help where it counts.** For fragments the model had already
studied, cutting at the seams did improve predictions. For fragments withheld
entirely, every version we built performed at chance — the same as guessing. The
improvement did not transfer.

**And most of the apparent benefit was not about the middle at all.** When we gave
a model nothing but the identity of the two inherited pieces — no middle sequence
whatsoever — it performed almost as well as the full method. So what looked like
insight into the creative, random middle was mostly bookkeeping about which
inherited parts the body happened to pick.

**A published neural network did no better.** The obvious objection is that our
methods were too simple. So we took a published deep-learning model built by
other researchers, retrained it on exactly the same data, and tested it the same
way. On its own practice material it scored 92%. On fragments it had never seen,
it landed at chance, like everything else. The problem is not that our tools were
too blunt.

## The check that changed one of our own numbers

Late in the work we asked a skeptical question about our own strongest result.

We identify the inherited front and back by looking for letters that barely vary
across thousands of receptors. But that rule stops a little early. Some inherited
letters end up counted as part of the "random" middle — and those letters are
perfectly predictable from which inherited piece was used. A model could look
clever simply by recognizing them, without understanding the middle at all.

So we deleted the letters at each edge of the middle and ran everything again.
The advantage dropped by about 60%. It did not vanish — real structure remains in
the interior, and it survives every control we could think of — but it is roughly
a third of what we first measured.

We report both numbers in the paper. The alternative was letting a reviewer, or a
later researcher, find it for us.

## Two side findings that may matter more than the main one

**A trap in public data.** For one of the three most-studied fragments, a single
laboratory contributed 79.5% of all the data. That lab's receptors carry an
unusual chemical signature — a particular amino acid appears in 10.8% of their
sequences against roughly 0.5% everywhere else. Anyone analyzing that fragment
without checking sources would "discover" a convincing biological signature that
is really one lab's experimental method. We nearly did exactly that, and we say so
in the paper.

**Counting the wrong thing.** The main public database lists the same receptor
once for every study that reported it. One sequence appears in 1,077 entries but
represents a single receptor. Count entries instead of receptors and you inflate
some patterns a thousandfold. An earlier version of our own analysis did this. We
found it, fixed it, and reported it.

## What we think it means

The receptor side of this problem is not the bottleneck.

There is a useful distinction hiding in our results: understanding *how receptors
are built* is not the same as understanding *what they recognize*. We got better
at the first and it bought us nothing on the second. Nothing in a description of
a receptor tells a model how a brand-new fragment should be matched against it.

So the effort probably needs to move to the other side of the lock and key —
better representations of the fragment and the molecule presenting it, rather than
ever more refined descriptions of the receptor. That is a redirection, and a
negative result is a legitimate way to argue for one.

## What this does not show

It does not show that the language idea is wrong in principle. It shows that this
particular implementation, at the scale of data currently available publicly, does
not transfer to new fragments. Richer approaches remain untested, including two
published methods that specifically claim to handle unseen fragments and that we
did not evaluate — we say so plainly in the paper rather than implying we covered
everything.

Our data also cover only one of the receptor's two arms, and our negative examples
are constructed by pairing receptors with fragments they probably don't bind,
rather than by experiment. Both limits are stated in the paper.

## Where the work lives

Everything is public: `https://github.com/etejkowski/grammar-of-immunity`. The
main analyses use nothing but the Python standard library, run in about six
minutes on a laptop, and produce identical numbers every time. Anyone can check
any figure in this summary against the code that produced it.
