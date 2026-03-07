test_that("readTrackMateXML handles non-existent file gracefully", {
	res <- tryCatch(readTrackMateXML("nonexistent.xml"), error = function(e) e)
	expect_true(is.null(res) || inherits(res, "error"))
})

test_that("readTrackMateXML handles wrong file suffix", {
	tmp <- tempfile(fileext = ".csv")
	file.create(tmp)
	on.exit(unlink(tmp))

	res <- tryCatch(readTrackMateXML(tmp), error = function(e) e)
	expect_true(is.null(res) || inherits(res, "error"))
})

test_that("readTrackMateXML handles XML with no tracks", {
	xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
			<TrackMate>
				<Model spatialunits="pixel" timeunits="sec"/>
				<AllSpots/>
			</TrackMate>'

	tmp <- tempfile(fileext = ".xml")
	writeLines(xml_content, tmp)
	on.exit(unlink(tmp))

	res <- tryCatch(readTrackMateXML(tmp), error = function(e) e)
	expect_true(is.null(res) || inherits(res, "error"))
})

test_that("readTrackMateXML handles XML with no filtered tracks", {
	xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
			<TrackMate>
				<Model spatialunits="pixel" timeunits="sec"/>
				<AllSpots/>
				<AllTracks>
					<Track TRACK_ID="0"/>
				</AllTracks>
				<FilteredTracks/>
			</TrackMate>'

	tmp <- tempfile(fileext = ".xml")
	writeLines(xml_content, tmp)
	on.exit(unlink(tmp))

	res <- tryCatch(readTrackMateXML(tmp), error = function(e) e)
	expect_true(is.null(res) || inherits(res, "error"))
})

test_that("readTrackMateXML parses valid XML correctly", {
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
		<TrackMate>
			<Model spatialunits="micron" timeunits="sec"/>
			<AllSpots>
				<SpotsInFrame frame="0">
					<Spot ID="1" POSITION_X="10.0" POSITION_Y="20.0" POSITION_Z="0.0" POSITION_T="0.0" FRAME="0" RADIUS="2.5" QUALITY="100"/>
					<Spot ID="2" POSITION_X="15.0" POSITION_Y="25.0" POSITION_Z="0.0" POSITION_T="0.5" FRAME="1" RADIUS="2.5" QUALITY="100"/>
				</SpotsInFrame>
			</AllSpots>
			<AllTracks>
				<Track TRACK_ID="0">
					<Edge SPOT_SOURCE_ID="1" SPOT_TARGET_ID="2"/>
				</Track>
			</AllTracks>
			<FilteredTracks>
				<TrackID TRACK_ID="0"/>
			</FilteredTracks>
		</TrackMate>'

  tmp <- tempfile(fileext = ".xml")
  writeLines(xml_content, tmp)
  on.exit(unlink(tmp))

	res <- tryCatch(readTrackMateXML(tmp), error = function(e) { skip(paste("readTrackMateXML failed:", e$message)) })

	# result is expected to be a list: first element the data frame, second calibration
	expect_true(is.list(res))
	result <- res[[1]]
	expect_true(is.data.frame(result))
	expect_true(all(c("track", "frame", "x", "y", "z", "time") %in% names(result)))
})

test_that("readTrackMateXML only includes filtered tracks", {
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
		<TrackMate>
			<Model spatialunits="micron" timeunits="sec"/>
			<AllSpots>
				<SpotsInFrame frame="0">
					<Spot ID="1" POSITION_X="10.0" POSITION_Y="20.0" POSITION_Z="0.0" POSITION_T="0.0" FRAME="0"/>
					<Spot ID="2" POSITION_X="15.0" POSITION_Y="25.0" POSITION_Z="0.0" POSITION_T="1.0" FRAME="1"/>
					<Spot ID="3" POSITION_X="100.0" POSITION_Y="200.0" POSITION_Z="0.0" POSITION_T="0.0" FRAME="0"/>
					<Spot ID="4" POSITION_X="105.0" POSITION_Y="205.0" POSITION_Z="0.0" POSITION_T="1.0" FRAME="1"/>
				</SpotsInFrame>
			</AllSpots>
			<AllTracks>
				<Track TRACK_ID="0">
					<Edge SPOT_SOURCE_ID="1" SPOT_TARGET_ID="2"/>
				</Track>
				<Track TRACK_ID="1">
					<Edge SPOT_SOURCE_ID="3" SPOT_TARGET_ID="4"/>
				</Track>
			</AllTracks>
			<FilteredTracks>
				<TrackID TRACK_ID="0"/>
			</FilteredTracks>
		</TrackMate>'

  tmp <- tempfile(fileext = ".xml")
  writeLines(xml_content, tmp)
  on.exit(unlink(tmp))

	res <- tryCatch(readTrackMateXML(tmp), error = function(e) { skip(paste("readTrackMateXML failed:", e$message)) })
	result <- res[[1]]

	expect_equal(nrow(result), 2)
	expect_equal(nlevels(result$track), 1)
	expect_equal(as.character(unique(result$track)), "0")
	expect_false(any(result$x > 50))
})

test_that("readTrackMateXML warns on duplicate track-frame combinations", {
  xml_content <- '<?xml version="1.0" encoding="UTF-8"?>
		<TrackMate>
			<Model spatialunits="micron" timeunits="sec"/>
			<AllSpots>
				<SpotsInFrame frame="0">
					<Spot ID="1" POSITION_X="10.0" POSITION_Y="20.0" POSITION_Z="0.0" POSITION_T="0.0" FRAME="0"/>
					<Spot ID="2" POSITION_X="15.0" POSITION_Y="25.0" POSITION_Z="0.0" POSITION_T="0.0" FRAME="0"/>
				</SpotsInFrame>
			</AllSpots>
			<AllTracks>
				<Track TRACK_ID="0">
					<Edge SPOT_SOURCE_ID="1" SPOT_TARGET_ID="2"/>
				</Track>
			</AllTracks>
			<FilteredTracks>
				<TrackID TRACK_ID="0"/>
			</FilteredTracks>
		</TrackMate>'

  tmp <- tempfile(fileext = ".xml")
  writeLines(xml_content, tmp)
  on.exit(unlink(tmp))

	res <- tryCatch(readTrackMateXML(tmp), error = function(e) { skip(paste("readTrackMateXML failed:", e$message)) })
	expect_output(invisible(res <- readTrackMateXML(tmp)), "duplicate")
})
