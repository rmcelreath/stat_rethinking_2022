<img src="title.gif" width="500">

Statistical Rethinking (2022 Edition)
===============

Instructor: Richard McElreath

Lectures: Uploaded and pre-recorded, two per week

Discussion: Online, Fridays 3pm-4pm Central European Time

# Purpose

This course teaches data analysis, but it focuses on scientific models first. The unfortunate truth about data is that nothing much can be done with it, until we say what caused it. We will prioritize conceptual, causal models and precise questions about those models. We will use Bayesian data analysis to connect scientific models to evidence. And we will learn powerful computational tools for coping with high-dimension, imperfect data of the kind that biologists and social scientists face.

# Format

Online, flipped instruction. The lectures are pre-recorded. We'll meet online once a week for an hour to work through the solutions to the assigned problems.

We'll use the 2nd edition of my book, <[Statistical Rethinking](https://xcelab.net/rm/statistical-rethinking/)>. I'll provide a PDF of the book to enrolled students.

Registration: Please sign up via <[COURSE IS FULL SORRY]>. I've also set aside 100 audit tickets at the same link, for people who want to participate, but who don't need graded work and course credit.

# Calendar & Topical Outline

There are 10 weeks of instruction. Links to lecture recordings will appear in this table. Weekly problem sets are assigned on Fridays and due the next Friday, when we discuss the solutions in the weekly online meeting.

[//]: # (11 Feb SPP conflict , 25 Feb Winter Break conflict )

| Week ## | Meeting date | Reading | Lectures |
| ------- | -------------- | ------------- | ---------------------- |
| Week 01 | 07 January  | Chapters 1, 2 and 3 | [1] <[The Golem of Prague](https://youtu.be/cclUd_HoRlo)> <[(Slides)](https://speakerdeck.com/rmcelreath/statistical-rethinking-2022-lecture-01)> <br> [2] Models & Bayesian Updating 
| Week 02 | 14 January | Chapter 4 | [3] Basic Regression <br> [4] Not-so-basic Regression
| Week 03 | 21 January | Chapters 5 and 6 |  [5] Confounding <br> [6] Even Worse Confounding
| Week 04 | 28 January | Chapters 7 and 8 | [7] Overfitting <br> [8] Interactions
| Week 05 | 04 February | Chapters 9, 10 and 11 | [9] Markov chain Monte Carlo <br> [10] Binomial GLMs
| Week 06 | 11 February | Chapters 11 and 12 | [11] Poisson GLMs <br> [12] Ordered Categories
| Week 07 | 18 February | Chapter 13 | [13] Multilevel Models <br> [14] Multi-Multilevel Models
| Week 08 | 25 February | Chapter 14 | [15] Varying Slopes <br> [16] Gaussian Processes
| Week 09 | 04 March | Chapter 15 | [17] Measurement Error <br> [18] Missing Data
| Week 10 | 11 March | Chapters 16 and 17 | [19] Beyond GLMs: State-space Models, ODEs <br> [20] Horoscopes


# Coding

This course involves a lot of scripting. Students can engage with the material using either the original R code examples or one of several conversions to other computing environments. The conversions are not always exact, but they are rather complete. Each option is listed below.

## Original R Flavor

For those who want to use the original R code examples in the print book, you need to install the `rethinking` R package. The code is all on github <https://github.com/rmcelreath/rethinking/> and there are additional details about the package there, including information about using the more-up-to-date `cmdstanr` instead of `rstan` as the underlying MCMC engine.

## R + Tidyverse + ggplot2 + brms

The <[Tidyverse/brms](https://bookdown.org/content/4857/)> conversion is very high quality and complete through Chapter 14.

## Python and PyMC3

The <[Python/PyMC3](https://github.com/pymc-devs/resources/tree/master/Rethinking_2)> conversion is quite complete.

## Julia and Turing

The <[Julia/Turing](https://github.com/StatisticalRethinkingJulia)> conversion is not as complete, but is growing fast and presents the Rethinking examples in multiple Julia engines, including the great <[TuringLang](https://github.com/StatisticalRethinkingJulia/TuringModels.jl)>.

## Other

The are several other conversions. See the full list at <https://xcelab.net/rm/statistical-rethinking/>.

# Homework and solutions

I will also post problem sets and solutions. Check the folders at the top of the repository.




