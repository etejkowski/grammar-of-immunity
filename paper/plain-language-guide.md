# The Grammar of Immunity — Plain Language Guide

A non-technical walkthrough of what this study asked, what we did, what the six
figures show, and what we actually found.

---

## The problem in one paragraph

Your immune system has to recognize germs it has never met — including ones that
don't exist yet. It solves this by manufacturing billions of different detectors
called **T-cell receptors**, each one shaped to grab a specific fragment of a
germ. If we could look at a receptor and predict *which* germ it grabs, we could
read someone's immune history from a blood sample, design vaccines faster, and
tell which cancer patients will respond to which treatment. Nobody can do this
reliably. Computer models get roughly nowhere when the target germ is one they
weren't trained on — which is exactly the situation in a new pandemic.

## The idea we tested

Receptors are built the way words are built.

Your body assembles each receptor from a small inventory of interchangeable
parts: a beginning piece, an end piece, and a randomly generated middle that
glues them together. That's very close to how words work in a language — a
prefix, a root, a suffix. English does this too: *un-* + *happy* + *-ness*.

A Nobel Prize winner named Niels Jerne noticed this in 1984 and titled his Nobel
lecture "The Generative Grammar of the Immune System." He proposed the idea but
never worked it out. In 2024 a team of linguists and immunologists revisited it
and argued that today's AI models fail at this problem because they chop receptor
sequences into meaningless chunks, ignoring the natural seams where the parts
join. They said someone should build a model that respects those seams. They
didn't build it.

**This study built it and tested it.**

The specific bet: if you cut receptor sequences at their natural joints — the way
you'd split *unhappiness* into *un-happy-ness* rather than into *unh-app-ine-ss*
— the model should learn something general about how receptors work, and
therefore handle germs it has never seen.

---

## What the pictures show

### Figure 1 — What "cutting at the seams" means

Five real receptor sequences, each shown as a row of letters (letters are amino
acids, the building blocks of proteins). Each row is colored in three parts:

- **Blue** — the beginning piece, inherited from your genes. Nearly identical in
  everyone.
- **Red** — the random middle. This is where each person's receptors differ, and
  it's the part that mostly decides what the receptor grabs.
- **Green** — the end piece, also inherited.

The point: what looks like one long meaningless string is really three
meaningfully different regions. Most computer models never make this distinction.
This figure is what our method does before anything else happens.

### Figure 2 — Yes, there is real structure in the random middle

We asked a simple question: is the "random" middle actually random?

Test: take a real middle region and shuffle its letters, like scrambling
*STRESSED* into *DESSERTS*. Same letters, same length, different order. Then see
if a model can tell the real one from the scrambled one. If it can, the order
carries information — the middles follow rules.

Four bars:

- **Gray, 0.50** — a deliberately useless model that only looks at length. It
  scores exactly 50%, i.e. pure coin-flip. This is a check that our test is
  honest: since shuffling doesn't change length, this model *must* fail, and it
  does, exactly.
- **Blue, 0.61** — a model that looks at letter order. Better than a coin flip.
  So the middles are not random.
- **Red, 0.72** — the same model, but told which inherited beginning and end
  pieces the receptor uses. Considerably better. Context matters.
- **Gray, 0.61** — the critical control. Same model, same complexity, but we
  *lied* to it about which pieces the receptor uses. Its advantage vanishes.

That last bar matters more than it looks. Bigger models often score better just
by being bigger. By feeding it scrambled information and watching the gain
disappear, we proved the improvement comes from real biology, not from the model
being fancier.

**Verdict: the linguistic structure is real and measurable.** So far, the idea is
working.

### Figure 3 — And it doesn't help at all where it counts

This is the important figure, and it's where the idea fails.

We tested on an official public competition dataset, using the competition's own
scoring rules, so we couldn't grade our own homework. The test contains 20 germ
fragments. Thirteen appear in the training data; **seven do not**.

Each method gets two bars: **solid** for the familiar targets, **striped** for
the unfamiliar ones.

Look at the solid bars — everything works reasonably. Look at the striped bars —
every single one collapses onto the coin-flip line. Including ours. Including
methods better than ours.

Cutting at the natural seams helps for germs the model has already studied. For
new germs, it buys nothing. And that's the only case anyone actually needs.

There's a subtle detail: two of the methods score *exactly* 50% on the unfamiliar
targets. That's not bad luck. Those methods work by looking up the most similar
receptor they've seen before, so when a germ is genuinely new they have nothing
to look up and simply shrug. They're not tuned badly — they're structurally
incapable of the task.

The rightmost pair of bars is the strongest version of the objection we could
build against ourselves. It's a published neural network (NetTCR-2.2) that reads
all six loops of both receptor arms, retrained by us on exactly the same training
data everything else here saw. If our failure were just a matter of using a
simple method, this is where it would show. It doesn't: the network matches our
crudest approach on familiar germs and lands at the coin-flip line on unfamiliar
ones. On its own practice questions — new receptors for germs it studied — it
scores 92%. It learned the material; it just can't transfer it to a germ it has
never seen.

### Figure 4 — Our one positive result turned out to be about data volume

We did find one real improvement: cutting at the seams beat naive chopping, for
familiar germs. Figure 4 asks whether that holds up.

Four bars, four ways of building the training set. The first two (red, positive)
used the full training data. The last two (gray, essentially zero) used about
*half* as much data.

Our first interpretation was that the harder test conditions killed the
advantage. The fourth bar exists to check that, and it disproves it: it has the
same reduced amount of data but easy conditions restored, and the advantage is
still gone.

So the improvement wasn't really about our clever cutting. It was about having
enough data. With plenty, the method looks better; with half, it doesn't. That's
a very different claim, and an honest paper has to make the smaller one.

There's a hopeful reading too: an effect that appears between "half" and "full"
might keep growing with far more data than we had.

### Figure 5 — A trap that catches a lot of published research

Early on we found what looked like a strong discovery: receptors targeting the
Epstein-Barr virus contained an unusual chemical building block (cysteine) far
more often than receptors targeting other viruses. It looked like a signature.

**Left panel:** where the data came from. For flu, no single laboratory
contributed more than about a third. For Epstein-Barr, **one laboratory
contributed 80%**.

**Right panel:** the "signature." Nearly 11% for Epstein-Barr versus under 1%
elsewhere.

Those two panels together explain the finding away. We weren't detecting the
virus. We were detecting one laboratory's equipment, methods, and habits. In
public databases of this kind, "which virus" and "which lab" are often the same
question, and a real-looking discovery can be nothing but the lab's fingerprint.

We had already written that finding up as a headline result. It's now in the
paper as a warning instead.

### Figure 6 — How we counted wrong, and the fix

**Left panel (note the sideways-stretched scale):** one particular receptor
sequence appears **1,077 times** in the database — but those are 1,077 database
entries describing what is really *one* receptor, reported over and over by
different studies. Counting entries instead of distinct receptors overstates how
common it is by a thousandfold.

We made exactly this mistake. A number in an early draft said "1,147 receptors"
when it should have said "database entries, mostly duplicates of one receptor."
We caught it only when we tried to draw this picture and the numbers came out
absurdly small.

**Right panel:** the same finding, counted properly. Flu-targeting receptors
really do share a distinctive pattern — 41% use the same inherited beginning
piece, and 13% carry a particular two-letter sequence in the middle, versus a few
percent for the other viruses.

Here's the thing, though: **this was published in 1998.** Our method rediscovered
a 28-year-old finding. That's genuinely good news — it shows the method sees real
biology — but it isn't a discovery, and we say so plainly.

### Figure 7 — The one bit of good news

If our clever cutting only works when there's lots of data, the obvious question
is: what happens as you add more? Does it help up to a point and then stop, or is
it still improving when you run out of data?

We trained the same models on six different amounts of data, from about 7,000
examples up to about 68,000, and plotted the result.

**Left panel:** the red line (our method) climbs steadily as data increases. The
blue line (naive chopping) is flat — it does not improve at all. More data simply
doesn't help the naive approach, while it keeps helping ours.

**Right panel:** the size of our advantage, plotted against how much data we
trained on. It rises the whole way and shows no sign of levelling off at the
largest amount we could get our hands on.

That's genuinely encouraging: the limit we ran into isn't a flaw in the method,
it's just that the world's public collection of this data is small. The method is
hungry, not broken.

The sobering part is in the purple line, the one for brand-new germs. It's also
rising — but so gently that if you follow the trend, you'd need roughly
**2,400 times** more data than exists today to get to a useful level. That's not
a matter of waiting a few years. It's the strongest argument in the whole study
for attacking the problem from the other direction.

---

## What we concluded

1. **The linguistic structure is real.** Receptor middles follow rules, and those
   rules depend on which inherited pieces they're attached to. This is solidly
   demonstrated, with the controls to back it.

2. **It doesn't transfer to new germs.** Nothing we or anyone else tried performs
   better than a coin flip on unfamiliar targets. The gap does close very
   slightly as we add data — but on that trend it would take thousands of times
   more data than exists to matter.

3. **The method is data-hungry, not capped.** Where it does help, it keeps
   helping as data grows, while the naive approach flatlines. That's a point in
   favor of the underlying idea.

4. **What little advantage exists is mostly explained by something simpler.**
   Just knowing which inherited pieces a receptor uses — no clever analysis of
   the random middle at all — performs about as well as our full method. The
   creative middle, which is the whole point of the linguistic idea, adds very
   little.

5. **We think the problem is on the other side.** Everything here describes the
   *receptor*. But nothing about a receptor tells you how a brand-new germ
   fragment maps onto the receptors that will catch it. We were carefully
   optimizing the half of the problem that isn't broken. The better question is
   what makes a germ fragment recognizable in the first place.

## Why publishing this is worth doing

It's a negative result, and negative results with proper controls are
scarce and useful. Specifically, it saves other people from spending a year on
the same reasonable-sounding idea, and it carries three warnings that apply well
beyond this project:

- Count distinct things, not database rows.
- Check whether your discovery is just one laboratory.
- Test at more than one data size, or you can't tell a better method from a
  hungrier one — and if you do sweep sizes, you learn which one you have.

Every claim here comes with the code that produces it, runs on a laptop in a few
minutes, and needs no special software.

## An honest note on how this went

Three separate times, a number that favored our idea turned out to be wrong when
checked more carefully — once from re-running with a different random starting
point, once from drawing a picture of it, once from adding a control. Each time
the corrected number was less exciting.

That's the normal shape of doing this properly. The results that survive that
kind of checking are the ones worth putting your name on.
