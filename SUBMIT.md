# Editorial Manager submission walkthrough

Everything below was checked against the journal's Guide for Authors, not
assumed. Where the Guide is silent, that is stated rather than guessed.

## Before you start

**This journal is not double-anonymized.** The Guide names no anonymized or
blinded review model anywhere, and it requires an "Essential title page
information" section *inside* the article. So author names and affiliations
stay in the manuscript and no separate anonymized file is needed. The
attachment's anonymization instructions do not apply here.

No filename in the upload set contains an author name.

## Article Type

The live drop-down lists 21 types, more than the Guide describes, and splits
the Guide's single "Full-length/ Research Article" into two separate entries:
**Full Length Article** and **Research Paper**. Either routes a complete report
of original research identically.

The last entry is **VSI: Autonomous Agents for Radiology Research** — the
special issue this manuscript was restructured for. Choose it if the call's
scope covers non-imaging biomedical work; the manuscript contains no radiology
content, deliberately, so if the call is strictly imaging-scoped take **Full
Length Article** instead and avoid an out-of-scope desk rejection.

## Upload these files

The manuscript goes up **first and alone** — Editorial Manager parses it for
title, abstract and references before the screen refreshes.

| Order | Item Type | File |
|---|---|---|
| 1 | Manuscript | `hf-crosstalk-manuscript.docx` |
| 2 | Title Page (with Author Details) | `Title Page.docx` |
| 3 | Declaration of Interest Statement | `Declaration of Interests.docx` |
| 4 | Highlights (for review) | `Highlights.docx` |
| 5 | *additional* | `Ethics Approval.docx` |
| 6 | *additional* | `CRediT Author Statement.docx` |
| 7 | *additional* | `figures_vector.zip` |

Every file is `.docx` except the figure archive. Editorial Manager rejects PDFs
as source files, and nothing in this set offers one.

**Graphical Abstract** is listed but optional. There isn't one, and the item can
be left empty.

### The declaration form

`Declaration of Interests.docx` is the journal's **own downloaded form**, filled
in — not a statement I wrote. The first of its three options is ticked ("no
known competing financial interests or personal relationships"), which is the
accurate one. The second option, for authors holding a journal role, is left
unticked because you confirmed no such relationship exists, and its placeholder
fields are cleared so no stray "Click here to enter your name" text survives.

Supersedes the earlier `Declaration of Interest.docx` I generated before the
form appeared. Upload the form, not that file.

### Title page and the manuscript

Editorial Manager wants the title page as its own item, so `Title Page.docx`
carries the title, all three authors with their affiliation letters, and the
corresponding author's full postal address and email.

Author details also remain inside the manuscript. That is correct here: this
journal is **not** double-anonymized, and its Guide requires an "Essential title
page information" section within the article. The duplication is intended, not
an oversight.

### Ethics

The screen asks for ethics approval or an explicit statement of
non-applicability. `Ethics Approval.docx` opens with the exact sentence they
specify — "Ethics Approval: Not applicable." — then states why: the study
re-analyzes only de-identified public datasets, with no human participants,
animals, or clinical trial.

### The figure archive

The archive unpacks on upload, so the eight PDFs land with no Item Type and show
in **red**. Use "Change Item Type of all" at the top of the list, set every
`.pdf` to Figure, and they clear in one action.

### About the figure archive

Editorial Manager unpacks archives on upload, so the eight PDFs will land in
the file list with no Item Type and appear in red. Use **"Change Item Type of
all"** at the top of the list, set every `.pdf` to **Figure**, and they clear
in one action. The filenames carry the figure numbers, so the order is
unambiguous.

The figures are vector PDFs with embedded fonts rather than raster images:
resolution-independent at any print size, which is what the Guide's artwork
table asks for. Each figure is a separate file, as required. Captions live in
the manuscript and are not attached to the figures — also as required, and
verified: no caption's title text appears inside its own figure.

## Field values for the later screens

**Full Title** (copy exactly):

> A multi-agent skill pipeline finds no tumor-selective heart-failure
> ligand–receptor axis across 30 tumor types under rank-based source criteria

**Abstract:** 239 words, structured as Background / Methods / Results /
Conclusions. If the parser mangles the structure, paste it from the
manuscript.

**Keywords** (the field wants them semicolon-separated):

> Cardio-oncology; heart failure; ligand-receptor signaling; multi-agent
> systems; computational reproducibility; pan-cancer transcriptomics

**Authors:** three, in order — Nicolas Avelar (corresponding, `a,b`), Ruidong
Zhang (`b`), Yi Pan (`b`). Affiliation `a` is North Oconee High School, Bogart,
GA; `b` is the School of Computing, University of Georgia, 200 D.W. Brooks
Drive, Athens, GA 30602. Editorial Manager will flag imported authors in red
until each has an email and affiliation.

**Region of Origin:** United States.

**Funding Information:** none. There is no grant to enter, and the Funding
section says so.

**Request Waiver:** the article-processing charge is waived for submissions
made before 31 December 2026, so there should be nothing to pay. If the fee
screen still appears, tick the waiver box and give the reason plainly — that
the work was self-funded with no institutional or grant support.

**Suggest / Oppose Reviewers:** optional. If you suggest anyone, the Guide
bars nominees from the same institution as any author, which rules out anyone
at the University of Georgia.

## Two things only you can do

- **Article Type** must be selected before anything else, and it changes which
  later steps appear. Pick the one matching a full research paper rather than a
  short communication or review.
- **Approve the built PDF.** The submission is not submitted until you open it
  from *Submissions Waiting for Approval by Author* and click Approve. Check
  the figures rendered and the Greek characters and superscripts survived; if
  anything is wrong, use Edit Submission rather than starting over.

## If something is rejected at upload

Send me the exact wording of the error. A file-format or word-count complaint
is usually a counting convention rather than a real overrun, and the gap tells
me which convention the portal used.
