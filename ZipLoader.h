#pragma once

#include <filesystem>
#include <zip.h>


class ZipLoader
{
public:
	explicit ZipLoader( const std::string& pathToArchive );
	~ZipLoader() noexcept;

	std::filesystem::path unarchive() const;

private:
	const std::filesystem::path pathToZip_;

    // TODO: make std::unique_ptr<zip>()
	zip* zipFile_ = nullptr; // zip file descriptor

};
