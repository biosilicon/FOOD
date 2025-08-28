import requests
from bs4 import BeautifulSoup

def fetch_td_texts(url):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, 'html.parser')

        divs = soup.find_all('div', {'data-v-ad57709a': True, 'class': 'column', 'data-testid': 'section-col'})

        td_texts = []
        for div in divs:
            tds = div.find_all('td')
            td_texts.extend(td.get_text(strip=True) for td in tds)

        return td_texts

    except requests.RequestException as e:
        print(f"Error fetching data from {url}: {e}")
        return []

if __name__ == "__main__":
    url = "https://nanoporetech.com/oo/organisms_seq"
    td_texts = fetch_td_texts(url)
    count = 1
    for text in td_texts:
        if count % 2 == 0:
            print(text)
        count += 1
        