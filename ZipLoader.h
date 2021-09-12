#pragma once

#include <filesystem>
#include <zip.h>


class ZipLoader
{
public:
	explicit ZipLoader( const std::string& pathToArchive );
	~ZipLoader() noexcept;

    zip* getZipFile() const;

	bool loadWaves();

private:
	const std::filesystem::path path_;

    // TODO: make std::unique_ptr<zip>()
	zip* zipFile_ = nullptr; // zip file descriptor

};
