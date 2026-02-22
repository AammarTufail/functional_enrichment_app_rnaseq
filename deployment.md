# Deployment Guide — Hostinger KVM VPS with Custom Domain

Complete step-by-step instructions to deploy the Functional Enrichment & GSEA App on a **Hostinger KVM VPS** with a custom domain, SSL, and Nginx reverse proxy.

---

## Prerequisites

- A **Hostinger KVM VPS** (minimum: 2 vCPU, 4 GB RAM, Ubuntu 22.04/24.04)
- A **domain name** pointed to Hostinger (e.g., `enrichment.yourdomain.com`)
- SSH access to the VPS
- Basic terminal knowledge

---

## Step 1 — Initial Server Setup

### 1.1 SSH into your VPS

```bash
ssh root@YOUR_VPS_IP
```

### 1.2 Update the system

```bash
apt update && apt upgrade -y
```

### 1.3 Create a non-root user (recommended)

```bash
adduser appuser
usermod -aG sudo appuser
su - appuser
```

### 1.4 Install essential packages

```bash
sudo apt install -y git curl wget nginx certbot python3-certbot-nginx ufw
```

---

## Step 2 — Configure Firewall (UFW)

```bash
sudo ufw allow OpenSSH
sudo ufw allow 'Nginx Full'
sudo ufw enable
sudo ufw status
```

---

## Step 3 — Point Your Domain to the VPS

### 3.1 In Hostinger DNS Zone Editor

Go to **Hostinger hPanel → Domains → DNS / Nameservers → DNS Records** and add:

| Type | Name | Value | TTL |
|------|------|-------|-----|
| A | `enrichment` (or `@` for root) | `YOUR_VPS_IP` | 3600 |

> If using a subdomain like `enrichment.yourdomain.com`, set **Name** to `enrichment`.
> If deploying on the root domain, set **Name** to `@`.

### 3.2 Verify DNS propagation

```bash
# Run from your local machine or the VPS
dig enrichment.yourdomain.com +short
# Should return YOUR_VPS_IP
```

> DNS propagation can take up to 24–48 hours, but usually resolves within minutes on Hostinger.

---

## Step 4 — Install Miniconda

```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
conda init bash
source ~/.bashrc
```

---

## Step 5 — Clone and Set Up the App

### 5.1 Clone the repository

```bash
cd ~
git clone https://github.com/AammarTufail/functional_enrichment_app_rnaseq.git
cd functional_enrichment_app_rnaseq
```

### 5.2 Create the conda environment

```bash
conda env create -f environment.yml
```

### 5.3 Test the app

```bash
conda activate enrichment_app
streamlit run app.py --server.port 8501 --server.headless true
# Press Ctrl+C to stop after confirming it works
```

---

## Step 6 — Create a Systemd Service

This ensures the app starts automatically on boot and restarts on failure.

### 6.1 Create the service file

```bash
sudo nano /etc/systemd/system/enrichment-app.service
```

Paste the following (adjust paths if needed):

```ini
[Unit]
Description=Functional Enrichment & GSEA Streamlit App
After=network.target

[Service]
User=appuser
Group=appuser
WorkingDirectory=/home/appuser/functional_enrichment_app_rnaseq
ExecStart=/home/appuser/miniconda3/envs/enrichment_app/bin/streamlit run app.py \
    --server.port 8501 \
    --server.headless true \
    --server.address 127.0.0.1 \
    --browser.gatherUsageStats false \
    --server.maxUploadSize 200
Restart=always
RestartSec=5
Environment="PATH=/home/appuser/miniconda3/envs/enrichment_app/bin:/usr/bin"

[Install]
WantedBy=multi-user.target
```

> **Note:** If you used `root` instead of `appuser`, change `User`, `Group`, and paths accordingly.

### 6.2 Enable and start the service

```bash
sudo systemctl daemon-reload
sudo systemctl enable enrichment-app
sudo systemctl start enrichment-app
sudo systemctl status enrichment-app
```

### 6.3 Useful service commands

```bash
# View logs
sudo journalctl -u enrichment-app -f

# Restart after code changes
sudo systemctl restart enrichment-app

# Stop the app
sudo systemctl stop enrichment-app
```

---

## Step 7 — Configure Nginx Reverse Proxy

### 7.1 Create Nginx config

```bash
sudo nano /etc/nginx/sites-available/enrichment-app
```

Paste (replace `enrichment.yourdomain.com` with your actual domain):

```nginx
server {
    listen 80;
    server_name enrichment.yourdomain.com;

    client_max_body_size 200M;

    location / {
        proxy_pass http://127.0.0.1:8501;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_read_timeout 86400;
        proxy_send_timeout 86400;
    }

    location /_stcore/stream {
        proxy_pass http://127.0.0.1:8501/_stcore/stream;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_read_timeout 86400;
    }
}
```

> The `/_stcore/stream` block ensures Streamlit's WebSocket connection works properly.

### 7.2 Enable the site

```bash
sudo ln -s /etc/nginx/sites-available/enrichment-app /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl reload nginx
```

### 7.3 Test HTTP access

Open `http://enrichment.yourdomain.com` in your browser — you should see the app.

---

## Step 8 — Enable HTTPS with Let's Encrypt

### 8.1 Obtain SSL certificate

```bash
sudo certbot --nginx -d enrichment.yourdomain.com
```

Follow the prompts:
- Enter your email
- Agree to Terms of Service
- Choose to redirect HTTP→HTTPS (recommended)

### 8.2 Verify auto-renewal

```bash
sudo certbot renew --dry-run
```

Certbot sets up a systemd timer for automatic renewal. Verify:

```bash
sudo systemctl status certbot.timer
```

### 8.3 Test HTTPS

Open `https://enrichment.yourdomain.com` — the app should load with a valid SSL certificate.

---

## Step 9 — Streamlit Configuration (Optional)

Create a Streamlit config file for production settings:

```bash
mkdir -p ~/functional_enrichment_app_rnaseq/.streamlit
nano ~/functional_enrichment_app_rnaseq/.streamlit/config.toml
```

```toml
[server]
headless = true
port = 8501
address = "127.0.0.1"
maxUploadSize = 200
enableCORS = false
enableXsrfProtection = true

[browser]
gatherUsageStats = false

[theme]
primaryColor = "#e74c3c"
backgroundColor = "#0e1117"
secondaryBackgroundColor = "#1a1a2e"
textColor = "#fafafa"
```

Restart the service after changes:

```bash
sudo systemctl restart enrichment-app
```

---

## Step 10 — Updating the App

When you push updates to the repository:

```bash
cd ~/functional_enrichment_app_rnaseq
git pull origin main
sudo systemctl restart enrichment-app
```

If dependencies changed:

```bash
conda activate enrichment_app
conda env update -f environment.yml --prune
sudo systemctl restart enrichment-app
```

---

## Troubleshooting

### App won't start

```bash
# Check service logs
sudo journalctl -u enrichment-app -n 50

# Test manually
conda activate enrichment_app
streamlit run app.py --server.port 8501 --server.headless true
```

### 502 Bad Gateway

- Verify the app is running: `sudo systemctl status enrichment-app`
- Check if port 8501 is listening: `ss -tlnp | grep 8501`
- Check Nginx logs: `sudo tail -f /var/log/nginx/error.log`

### WebSocket errors in browser

Ensure the `/_stcore/stream` location block is in your Nginx config with the `Upgrade` and `Connection` headers.

### SSL certificate issues

```bash
sudo certbot renew --force-renewal
sudo systemctl reload nginx
```

### Out of memory

If the VPS runs out of RAM during GSEA analysis, add swap:

```bash
sudo fallocate -l 4G /swapfile
sudo chmod 600 /swapfile
sudo mkswap /swapfile
sudo swapon /swapfile
echo '/swapfile none swap sw 0 0' | sudo tee -a /etc/fstab
```

---

## Architecture Overview

```
Browser (HTTPS)
    │
    ▼
Nginx (port 443) ──SSL──► Let's Encrypt
    │
    ▼ proxy_pass
Streamlit (127.0.0.1:8501)
    │
    ▼ KEGG REST API
https://rest.kegg.jp
```

---

## Security Checklist

- [x] Non-root user for running the app
- [x] UFW firewall enabled (SSH + Nginx only)
- [x] HTTPS via Let's Encrypt with auto-renewal
- [x] Streamlit bound to `127.0.0.1` (not exposed directly)
- [x] XSRF protection enabled
- [ ] (Optional) Add HTTP Basic Auth via Nginx for private access
- [ ] (Optional) Set up fail2ban for SSH brute-force protection

### Optional: Add password protection

```bash
sudo apt install apache2-utils
sudo htpasswd -c /etc/nginx/.htpasswd youruser
```

Add to the Nginx `location /` block:

```nginx
auth_basic "Restricted Access";
auth_basic_user_file /etc/nginx/.htpasswd;
```

```bash
sudo systemctl reload nginx
```
