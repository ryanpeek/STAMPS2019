# nicolas cage + movies


# Load Libraries ----------------------------------------------------------

library(tidyverse)
library(scales)
library(glue)


# Load Data From imdb -----------------------------------------------------


# https://jmgirard.com/tag/tidyverse/
# https://minimaxir.com/2018/07/imdb-data-analysis/

# see here:
# https://www.imdb.com/interfaces/

# get data:
# df_all <- read_tsv('https://datasets.imdbws.com/title.basics.tsv.gz', na = "\\N", quote = '')
#load("data/imdb_title_basics.rda")

#df_ratings <- read_tsv('https://datasets.imdbws.com/title.ratings.tsv.gz', na = "\\N", quote = '')
#load("data/imdb_title_ratings.rda")

#df_cast <- read_tsv("https://datasets.imdbws.com/title.principals.tsv.gz", na="\\N", quote='') %>% filter(str_detect(category, "actor|actress"))
#load("data/imdb_title_cast.rda")

# df_names <- read_tsv("https://datasets.imdbws.com/name.basics.tsv.gz", na = "\\N", quote = '') 
# df_names <- df_names %>% 
#   filter(str_detect(primaryProfession, "actor|actress"))  %>%
#   select(nconst, primaryName, birthYear, knownForTitles) 
#load("data/imdb_actor_ress_names.rda")

# filter to movies by nicolas cage
# imdb_sel <- df_all %>% 
#   dplyr::filter(
#     titleType == "movie",  # exclude non-movies
#     isAdult == 0,          # exclude adult movies
#     runtimeMinutes <= 240, # exclude movies over 4 hours long
#     startYear >= 1918     # exclude movies more than 100 years old
#   ) %>% 
#   dplyr::select(tconst, primaryTitle, startYear, runtimeMinutes, genres) %>% 
#   dplyr::arrange(startYear, primaryTitle) %>% 
#   print()
# 
# # join with cast
# df_cage <- df_names %>% filter(str_detect(primaryName, "Nicolas Cage"))
# df_cage_movies <- inner_join(df_cast, df_cage, by=c("nconst"))
# 
# # join with ratings
# df_cage_movies_rated <- left_join(df_cage_movies, df_ratings)
# 
# # rm big stuff
# rm(df_all, df_cast, df_names, df_ratings, imdb_sel, df_cage, df_cage_movies)
# 
# # save nic cage
# save(df_cage_movies_rated, file = "data/nicolas_cage_movies_rated.rda")


# START HERE --------------------------------------------------------------
library(tidyverse)
library(scales)
library(glue)

load("data/nicolas_cage_movies_rated.rda")


# Visualize the distribution of runtimes 
(t1 <- glue("Mean runtime (min) distribution for Nicolas Cage was {round(mean(df_cage_movies_rated$runtimeMinutes, na.rm=T), 2)} for movies in IMDB (1918-2019)"))

(t2 <- glue("Runtime (min) distribution for Nicolas Cage was {round(sum(!is.na(df_cage_movies_rated$runtimeMinutes)), 2)} for movies in IMDB (1918-2019)"))

df_cage_movies_rated %>% 
  ggplot(aes(x = runtimeMinutes)) +
  geom_histogram(binwidth = 5, fill = "maroon", color = "black") +
  scale_x_continuous(breaks = seq(0, 240, 60)) +
  scale_y_continuous(
    labels = scales::unit_format(scale = 1e-3, suffix = "k")) +
  labs(
    x = "Runtime (minutes)", 
    y = "Count (of movies)", 
    title = t2
  )

# avg rating
plotly::ggplotly(
  ggplot(df_cage_movies_rated, aes(x = numVotes, y = averageRating)) +
    geom_bin2d() +
    scale_x_log10(labels = comma) +
    scale_y_continuous(breaks = 1:10) +
    scale_fill_viridis_c(labels = comma)
)  

# split out the genre
df_cage_movies_rated <- df_cage_movies_rated %>% separate(genres, into = c("genre1", "genre2"))

(p1 <- ggplot(df_cage_movies_rated %>% 
         filter(titleType == "movie", grepl("Action|Comedy", genre1)), 
       aes(x = averageRating, y = runtimeMinutes, fill=genre1, label=primaryTitle)) +
  geom_point(pch=21, size=4) +
  #scale_y_continuous(breaks = seq(0, 180, 60), labels = 0:3) +
  scale_x_continuous(breaks = 0:10) +
  geom_smooth(method = "glm", aes(fill=genre1)) +
  scale_fill_viridis_d(option = "inferno") +
  #scale_fill_viridis_c("Avg Rating", option = "inferno", labels = comma) +
  theme_minimal(base_family = "Roboto Condensed", base_size = 8) +
  labs(title = "Relationship between Nicolas Cage Movie Runtime and Average Mobie Ratings",
       subtitle = "Data from IMDb retrieved July 4th, 2018",
       y = "Runtime (minutes)",
       x = "Average User Rating",
       caption = "RP-ryanpeek.github.com",
       fill = "Genre"))

plotly::ggplotly(p1)

# are nick cage action movies higher rated than his comedy movies?

library(tidymodels)

# Confidence interval for a difference in means (using the non-formula interface giving both the response and explanatory variables in specify()):
df_mod <- df_cage_movies_rated %>% filter(titleType=="movie", grepl("Action|Comedy", genre1)) %>% mutate(genre1=as.factor(genre1))

df_mod %>% 
  specify(response = averageRating, explanatory = genre1) %>%
  generate(reps = 100, type = "bootstrap") %>%
  calculate(stat = "diff in means", order = c("Action", "Comedy"))


library(brms)

mod_eqvar <- brm(
  averageRating ~ genre1, 
  data = df_mod)

summary(mod_eqvar)

# unequal var per group
mod_robust <- brm(bf(averageRating ~ genre1 + runtimeMinutes, 
     sigma ~ genre1),
  family=student,
  data = df_mod, 
  cores=4)

plot(mod_robust)
plot(mod_robust, N = 2, ask = FALSE)
plot(marginal_effects(mod_robust), points = TRUE)

get_variables(mod_robust)

# get draws
mod_robust %>%
  spread_draws(b_Intercept, b_sigma_Intercept) %>%
  head(10)

library(modelr)
df_mod %>%
  data_grid(genre1, runtimeMinutes) %>%
  add_fitted_draws(mod_robust) %>%
  ggplot(aes(x = .value, y = genre1)) +
  stat_pointintervalh(.width = c(.66, .95))



# hypothesis testing

hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_genre1Comedy) = 0")
hypothesis(mod_robust, hyp)


hyp <- "exp(sigma_Intercept + sigma_genre1Comedy) > exp(sigma_Intercept)"
(hyp <- hypothesis(mod_robust, hyp))

plot(hyp, chars = NULL)


glm1 <- glm(data = df_mod, averageRating ~ genre1)
summary(glm1)


library(tidybayes)
