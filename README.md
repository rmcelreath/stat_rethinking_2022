<img src="title.gif" width="500">

Statistical Rethinking (2022 Edition)
===============

Instructor: Richard McElreath

Lectures: Uploaded <[Playlist](https://www.youtube.com/playlist?list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN)> and pre-recorded, two per week

Discussion: Online, Fridays 3pm-4pm Central European Time

# Purpose

This course teaches data analysis, but it focuses on scientific models first. The unfortunate truth about data is that nothing much can be done with it, until we say what caused it. We will prioritize conceptual, causal models and precise questions about those models. We will use Bayesian data analysis to connect scientific models to evidence. And we will learn powerful computational tools for coping with high-dimension, imperfect data of the kind that biologists and social scientists face.

# Format

Online, flipped instruction. The lectures are pre-recorded. We'll meet online once a week for an hour to work through the solutions to the assigned problems.

We'll use the 2nd edition of my book, <[Statistical Rethinking](https://xcelab.net/rm/statistical-rethinking/)>. I'll provide a PDF of the book to enrolled students.

Registration: Please sign up via <[COURSE IS FULL SORRY]>. I've also set aside 100 audit tickets at the same link, for people who want to participate, but who don't need graded work and course credit.

# Calendar & Topical Outline

There are 10 weeks of instruction. Links to lecture recordings will appear in this table. Weekly problem sets are assigned on Fridays and due the next Friday, when we discuss the solutions in the weekly online meeting.

Lecture playlist on Youtube: <[Statistical Rethinking 2022](https://www.youtube.com/playlist?list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN)>

[//]: # (11 Feb SPP conflict , 25 Feb Winter Break conflict )

| Week ## | Meeting date | Reading | Lectures |
| ------- | -------------- | ------------- | ---------------------- |
| Week 01 | 07 January  | Chapters 1, 2 and 3 | [1] <[The Golem of Prague](https://youtu.be/cclUd_HoRlo)> <[(Slides)](https://speakerdeck.com/rmcelreath/statistical-rethinking-2022-lecture-01)> <br> [2] <[Bayesian Inference](https://www.youtube.com/watch?v=guTdrfycW2Q&list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN&index=2)> <[(Slides)](https://speakerdeck.com/rmcelreath/statistical-rethinking-2022-lecture-02)> 
| Week 02 | 14 January | Chapters 4 and 5 | [3] <[Basic Regression](https://www.youtube.com/watch?v=zYYBtxHWE0A)> <[(Slides)](https://speakerdeck.com/rmcelreath/statistical-rethinking-2022-lecture-03)> <br> [4] <[Categories & Curves](https://youtu.be/QiHKdvAbYII)> <[(Slides)](https://speakerdeck.com/rmcelreath/statistical-rethinking-2022-lecture-04)>
| Week 03 | 21 January | Chapters 5 and 6 |  [5] <[Elemental Confounds](https://youtu.be/UpP-_mBvECI)> <[(Slides)](https://speakerdeck.com/rmcelreath/statistical-rethinking-2022-lecture-05)> <br> [6] <[Good & Bad Controls](https://www.youtube.com/watch?v=NSuTaeW6Orc&list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN&index=6)> <[(Slides)](https://speakerdeck.com/rmcelreath/statistical-rethinking-2022-lecture-06)>
| Week 04 | 28 January | Chapters 7, 8 and 9 | [7] <[Overfitting](https://www.youtube.com/watch?v=odGAAJDlgp8&list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN&index=7)> <[(Slides)](https://speakerdeck.com/rmcelreath/statistical-rethinking-2022-lecture-07)> <br> [8] <[Markov chain Monte Carlo](https://www.youtube.com/watch?v=Qqz5AJjyugM&list=PLDcUM9US4XdMROZ57-OIRtIK0aOynbgZN&index=8&pp=sAQB)> <[(Slides)](https://speakerdeck.com/rmcelreath/statistical-rethinking-2022-lecture-08)>
| Week 05 | 04 February | Chapters 10 and 11 | [9] Generalized Linear Models <br> [10] Binomial GLMs
| Week 06 | 11 February | Chapters 11 and 12 | [11] Poisson GLMs <br> [12] Ordered Categories
| Week 07 | 18 February | Chapter 13 | [13] Multilevel Models <br> [14] Multi-Multilevel Models
| Week 08 | 25 February | Chapter 14 | [15] Varying Slopes <br> [16] Gaussian Processes
| Week 09 | 04 March | Chapter 15 | [17] Measurement Error <br> [18] Missing Data
| Week 10 | 11 March | Chapters 16 and 17 | [19] Beyond GLMs: State-space Models, ODEs <br> [20] Horoscopes


# Coding

This course involves a lot of scripting. Students can engage with the material using either the original R code examples or one of several conversions to other computing environments. The conversions are not always exact, but they are rather complete. Each option is listed below. I also list conversions <[here](https://xcelab.net/rm/statistical-rethinking/)>.

## Original R Flavor

For those who want to use the original R code examples in the print book, you need to install the `rethinking` R package. The code is all on github <https://github.com/rmcelreath/rethinking/> and there are additional details about the package there, including information about using the more-up-to-date `cmdstanr` instead of `rstan` as the underlying MCMC engine.

## R + Tidyverse + ggplot2 + brms

The <[Tidyverse/brms](https://bookdown.org/content/4857/)> conversion is very high quality and complete through Chapter 14.

## Python: PyMC3 and NumPyro and more

The <[Python/PyMC3](https://github.com/pymc-devs/resources/tree/master/Rethinking_2)> conversion is quite complete. There are also at least two NumPyro conversions: <[NumPyro1](https://github.com/asuagar/statrethink-course-numpyro-2019)> <[NumPyro2](https://fehiepsi.github.io/rethinking-numpyro/)>. And there is this <[TensorFlow Probability](https://github.com/ksachdeva/rethinking-tensorflow-probability)>. 

## Julia and Turing

The <[Julia/Turing](https://github.com/StatisticalRethinkingJulia)> conversion is not as complete, but is growing fast and presents the Rethinking examples in multiple Julia engines, including the great <[TuringLang](https://github.com/StatisticalRethinkingJulia/TuringModels.jl)>.

## Other

The are several other conversions. See the full list at <https://xcelab.net/rm/statistical-rethinking/>.

# Homework and solutions

I will also post problem sets and solutions. Check the folders at the top of the repository.




