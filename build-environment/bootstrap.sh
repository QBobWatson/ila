#!/usr/bin/env bash

set -e

cd /base/build-environment

# Fix broken resolved configuration
rm /etc/systemd/resolved.conf
systemctl restart systemd-resolved

# Provision the build environment

read -d '' PACKAGES <<EOF || true
    python-minimal
    python3-minimal
    python-pip
    python3-pip
    libcairo2-dev
    pkg-config
    trang
    xsltproc
    ruby-dev
    scons
    python-mako
    coffeescript
    npm
    apache2
    python-poppler
    fontforge
    python-fontforge
    texlive-full
    libxml2-utils
    inkscape
EOF

PYTHON2_PACKAGES=pdfrw

read -d '' PYTHON3_PACKAGES <<EOF || true
    cssutils
    lxml
    bs4
    mako
    pyyaml
EOF

apt-get update
apt-get install -y $PACKAGES

pip2 install $PYTHON2_PACKAGES
pip2 install pycairo-1.15.4.tar.gz

pip3 install $PYTHON3_PACKAGES

# Patched inkscape
dpkg -i inkscape-patched_0.92.3-1_amd64.deb

gem install compass

# Apache
rm -rf /var/www/html
ln -s /home/vagrant/build /var/www/html
