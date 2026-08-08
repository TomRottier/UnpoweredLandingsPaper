#set page(numbering: "1", margin: 0.5in)
#set heading(numbering: "1.")
#set text(12pt, font: "New Computer Modern")
#set math.equation(numbering: "(1)")
#set par(justify: true)

#align(center, text(18pt, weight: "bold")[
  Flight dynamics and optimal control of unpowered avian landings])

#align(center, text(13pt)[
  Tom Rottier, Philippe Lavoie, Ben Parslew
])

#include "abstract.typ"

#include "introduction.typ"

#include "methods.typ"

= Results and Discussion
#include "results-dynamics.typ"
#include "results-opt-control.typ"
#include "results-constraints.typ"
#include "results-perching.typ"
#include "results-sensitivity.typ"

#include "conclusion.typ"


#pagebreak()

#bibliography("references.bib", style: "proceedings-of-the-royal-society-b.csl")
