import { useRelevantArticles } from '../../hooks/useRelevantArticles';
import { Books, LinkSimple } from '@phosphor-icons/react';
import { Accordion, Button } from 'react-bootstrap';
import { highlightWords } from '../../utils/utils';
import { useUrlParams } from '../../hooks/useUrlParams';

const RelevantArticles = () => {
  const { relevantArticles, relevantArticlesError } = useRelevantArticles();
  const {
    params: { terms },
  } = useUrlParams();

  if (relevantArticlesError || relevantArticles.length === 0) {
    return;
  }

  return (
    <div className={'mb-5'}>
      <h3 className={'d-flex justify-content-center align-items-center mb-4'}>
        <Books className={'text-secondary'} weight={'light'} />
        <div className={'vr mx-2'} />
        Relevant Articles
      </h3>
      <Accordion>
        {relevantArticles.map((article, index) => (
          <Accordion.Item key={index} eventKey={index.toString()}>
            <Accordion.Header>
              <div>
                <div className={'mb-2'}>{article.title}</div>
                <div className={'d-flex align-items-center small text-secondary'}>
                  {article.authors && article.authors.join(', ')}
                  <div className={'vr mx-2'} />
                  {article.country}
                </div>
              </div>
            </Accordion.Header>
            <Accordion.Body>
              <div>{highlightWords(article.abstract, terms)}</div>
              <div className={'d-flex justify-content-end mt-3'}>
                <Button href={article.url} target={'_blank'} variant={'outline-primary'}>
                  <div className={'d-flex align-items-center'}>
                    <LinkSimple weight={'light'} />
                    <div className={'vr mx-2'} />
                    Read on PubMed
                  </div>
                </Button>
              </div>
            </Accordion.Body>
          </Accordion.Item>
        ))}
      </Accordion>
    </div>
  );
};

export default RelevantArticles;
