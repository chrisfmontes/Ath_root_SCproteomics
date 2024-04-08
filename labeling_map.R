library(dplyr)

# original code
set.seed(1)
labeling_map <- data.frame(plate = rep(paste0("plate_",1:2), each = 384),
                           pos = rep(paste0(rep(LETTERS[1:16], each = 24),1:24), times = 2)) |>
  slice(-(rep((0:1*384)+1,each = 6)+0:5)) |>
  group_by(plate) |>
  mutate(PlexNumber = sample(paste0(rep("Plex_",378),1:54), size = 378, replace = FALSE)) |>
  ungroup() |> group_by(PlexNumber) |>
  mutate(TMTlabel = sample(x = c(5:18), size = 14, replace = FALSE))
  
write.csv(x = labeling_map, file = "TMTlabelMap.csv")

# Let's make position matrices for the OT-2
for (i in 5:18) {
  labeling_map |>
    ungroup() |>
    filter(plate == "plate_1" & TMTlabel == i) |>
    select(pos) |>
    rename(!!paste0("TMT",i):=pos) |> t() |>
    write.table(file = "TMT_label_well_mapping_plate1.csv", append = TRUE, col.names = FALSE, sep = ",")
}

### Now we mark the wells to pool samples after labeling and make position matrices for OT-2
for (i in 1:54) {
  labeling_map |>
    ungroup() |>
    filter(PlexNumber == paste0("Plex_",i)) |>
    select(pos) |>
    rename(!!paste0("PlexNumber",i):=pos) |> t() |>
    write.table(file = "TMT_pool_well_mapping.csv", append = TRUE, col.names = FALSE, sep = ",")
}

# Specify wells by plate
for (i in 1:54) {
  labeling_map |>
    ungroup() |>
    filter(plate == "plate_2" & PlexNumber == paste0("Plex_",i)) |>
    select(pos) |>
    rename(!!paste0("PlexNumber",i):=pos) |> t() |>
    write.table(file = "TMT_pool_well_mapping_plate2.csv", append = TRUE, col.names = FALSE, sep = ",")
}
coord = paste0(rep(LETTERS[1:8], each = 12),1:12)
write.table(x = t(coord[1:54]), file = "96-well_coords_54Plex.csv", row.names = FALSE, col.names = FALSE, sep = ",")

            