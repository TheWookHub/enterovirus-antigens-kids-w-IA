# function to rename elements in a named vector based on substring matching
# use fixed in str_detect to use literal matching rather than regex matching 
rename_by_substring <- function(named_vector, name_map) {
  new_names <- map_chr(names(named_vector), function(name) {
    matched_key <- names(name_map)[str_detect(name, fixed(names(name_map)))]
    name_map[matched_key]
  })
  names(named_vector) <- new_names
  named_vector
}

### usage 
# test_some_named_vector <- c(a = 3, b = 9)
# test_name_map <- c("a" = "apple", "b" = "banana")

# rename_by_substring(test_some_named_vector, test_name_map)